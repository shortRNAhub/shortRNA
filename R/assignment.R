#' assignReads
#'
#' Assigns each read to one or more of possible sources.
#'
#' @param sources A data.frame of possible sources, as produced by the findOverlap function.
#' @param rules A list of assignment rules; defaults to `defaultAssignRules()`.
#' @param nthreads The number of threads to use; default to detected cores -1 if the parallel package is installed.
#'
#' @return A data.frame of assigned sources.
#'
#' @export
assignReads <- function(sources, rules=defaultAssignRules(), nthreads=NULL){
  library(data.table)
  if(!is(sources, "data.table")) sources <- data.table(sources,stringsAsFactors=F)

  nn <- c("seq", "cigar", "chrPos", "strandRead", "strandFeature", "posInFeature", 
         "overlap", "percentOverlap", "transcript_id", "gene_id", "transcript_type")
  if(!all(colnames(sources)==nn)) stop(paste("`sources` should have exactly the following column names:", paste(nn,collapse=", ")))
  nn <- c(nn, "length", "validType", "status")
  
  sources$length <- sapply(sources$seq, FUN=nchar)
  
  # we check which feature types are valid with the overlaping reads
  sources$validType <- sapply(split(sources,1:nrow(sources),drop=F),rules=rules,FUN=.validateFeatureType)
  # in the case of primary piRNA, invalid still gets assigned to the precursor
  sources[which(sources$validType & sources$transcript_type=="piRNA_precursor"),"transcript_type"] <- "primary_piRNA"
  sources[which(!sources$validType & sources$transcript_type=="piRNA_precursor"),"validType"] <- TRUE
  
  DEBUG <- FALSE
  
  if(DEBUG){
    par <- NULL
  }else{
    par <- .checkPara(nthreads)
  }
  
  if(is.null(par)){
    sources <- as.data.frame(rbindlist(lapply(split(sources,sources$seq,drop=F), rules=rules, DEBUG=DEBUG, FUN=function(sources,rules,DEBUG){ 
      a <- try(assignRead(sources,rules),silent=T)
      if(is(a,"try-error")){
        if(DEBUG){
          print("Error trying to assign:")
          print(sources)
        }
        stop(a[[1]])
      }
      return(a)
    })))
    colnames(sources) <- nn
    
  }else{
    fns <- c("assignRead", ".validateFeatureType", ".aggSources", ".agCigar", ".disambiguate_miRNAs", ".disambiguate_tRNA", "defaultAssignRules")
    sources <- foreach(x=split(sources,sources$seq,drop=F), .export=fns) %dopar% {
      res <- as.data.frame(assignRead(x,rules))
      names(res) <- nn
      res
    }
    sources <- as.data.frame(rbindlist(sources))
    stopCluster(par)
  }
  row.names(sources) <- sources[,1]
  sources <- sources[,-1]
  sources$validType <- NULL
  
  .recastSources(sources)
  
}

#' assignRead
#'
#' Assigns one read to one or more of possible sources.
#'
#' @param sources A data.frame of possible sources for the read (i.e., a subset of the `sources` slot of a shortRNAexp object).
#' @param rules A list of assignment rules; defaults to `defaultAssignRules()`.
#'
#' @return A vector including the following variables: seq, cigar, chrPos, strandRead, strandFeature, posInFeature, overlap, percentOverlap, transcript_id, gene_id, transcript_type, status
#'
#' @export
assignRead <- function(sources, rules=defaultAssignRules()){
    nn <- c("seq", "cigar", "chrPos", "strandRead", "strandFeature", "posInFeature", 
         "overlap", "percentOverlap", "transcript_id", "gene_id", "transcript_type", "length", "validType")
    if(!all(nn %in% names(sources))) stop(paste("`sources` should have at least the following column names:", paste(nn,collapse=", ")))
    
    # if there are multiple mapping locations,
    # we first disambiguate each location before working on each
    nloc <- length(unique(sources$chrPos))
    if(nloc>1 & nloc<nrow(sources)){
        sources <- rbindlist(lapply(split(sources,sources$chrPos,drop=F),rules=rules, FUN=assignRead))
    }

    # we eliminate overlaps for which the sequence is incompatible with the feature:
    w <- which(sources$validType)
    if(length(w)==0){
        # there is no overlap for which the read is valid
        if(nrow(sources)>1) sources <- .aggSources(sources)
        sources$transcript_type <- NA_character_
        sources$status <- "invalid"
        return(sources)
    }else{
        # there are valid overlaps
        if(nloc==1){
            # if all overlaps are from the same chromosomal position,
            # we eliminate the invalid ones
            sources <- sources[w,,drop=F]
        }else{
            # otherwise we simply flag the invalid ones as NAs
            sources$transcript_type[-w] <- NA_character_
        }
    }
    
    if(nloc==1){
        # we're dealing with overlaps from a single genomic location
        
        # prioritize features on the same strand as the read
        if(nrow(sources)>1 & rules$sameStrand=="prioritize"){
            w <- which(sources$strandRead==sources$strandFeature)
            if(length(w)>0){
                sources <- sources[w,,drop=F]
            }
        }
    
        # we eliminate the low-priority types if there's something else
        if(!all(sources$transcript_type %in% rules$lowPriorityTypes)){
            sources <- sources[which(!(sources$transcript_type %in% rules$lowPriorityTypes)),]
        }
    }

    if(nrow(sources)>1 & rules$prioritizeKnown){
        # if in the rules, we prioritize known over unknown sources; in other words,
        # if a read maps to multiple location and some of these overlap with a feature,
        # we eliminate those alignments that do not overlap any feature
        known <- !is.na(sources$transcript_id) & !is.na(sources$gene_id)
        if(all(!known)){
            # there is no overlap with known features, so we can report this
            sources <- sources[1,,drop=F]
            if(nloc>1) sources$chrPos <- "ambiguous"
            sources$status <- "unknown"
            return(sources)
        }else{
            # there are some known features overlapping, so we eliminate the 
            # alignments without known features
            sources <- sources[which(known),,drop=F]
        }
    }
    
    # if there are still multiple locations, we consider location ambiguous:
    if(length(unique(sources$chrPos))>1) sources$chrPos <- "ambiguous"
    sources <- sources[!duplicated(sources),]
  
    # we check if the sequence could be a secondary piRNA
    # ###################################3
    # TO DO
    # ###################################3
  
    # if there is any (non-eliminated) alignment without feature, we report as ambiguous:
    if(any(is.na(sources$transcript_type))){
        sources <- .aggSources(sources)
        sources$transcript_type <- NA_character_
        sources$status <- "ambiguous"
        return(sources)
    }
 
    # if there are multiple source features, we proceed with the one(s) that have the largest overlap
    if(nrow(sources)>1){
        s2 <- sources[which(sources$overlap==max(sources$overlap) |
        (sources$transcript_type %in% rules$highPriorityTypes & sources$overlap==max(sources$overlap)-1) | 
        (sources$transcript_type=="miRNA" & grepl("-(3|5)p",sources$transcript_id))),,drop=F]
        if(nrow(s2)>0) sources <- s2
    }

    # if there is any source type withing the high-priority list, we assign the read to it
    if(any(sources$transcript_type %in% rules$highPriorityTypes)){
        sources <- sources[which(sources$transcript_type %in% rules$highPriorityTypes),,drop=F]
    }

    # if there's a unique source, we return it:
    if(nrow(sources)==1){
        # we get the more specific tRNA fragment type:
        if(sources$transcript_type=="tRNA") sources$transcript_type <- tRFtype(sources, rules)
        if((is.na(sources$transcript_id) && is.na(sources$gene_id))){
            sources$status <- "unknown"
        }else{
            if(!("status" %in% names(sources))) sources$status <- "unique"
        }
        return(sources)
    }

    # if all overlaps are of the same type, we check if we can find a common meta-feature
    if(all(sources$transcript_type == "miRNA")) return(.disambiguate_miRNAs(sources))
    if(all(sources$transcript_type %in% c("pseudo_tRNA","tRNA")))   return(.disambiguate_tRNA(sources, rules))    
    
    # if all the overlaps point to the same feature, we're done:
    if(length(unique(sources$gene_id))==1){
        return(.aggSources(sources, addStatus=ifelse(all(is.na(sources$transcript_id)),"unknown","unique")))
    }
    
    # eventually add, here, a step of EM...
    
    if(length(unique(sources$transcript_type))>1) sources$transcript_type <- "ambiguous"
    sources <- .aggSources(sources)
    sources$status <- "ambiguous"
    return(sources)
}

check_miRNA <- function(src, length=20:23, minOverlap=12){
    if(nrow(src)>1) src <- src[1,,drop=F]
    if(!(src$length %in% length)) return(FALSE)
    if(src$overlap < minOverlap) return(FALSE)
    return(TRUE)
}

.disambiguate_miRNAs <- function(sources){
    sources$transcript_id <- gsub("p(A|C|G|T)$","p",sources$transcript_id)
    nn <- gsub("-([0-9a-b]+)-([0-9]*)$","-\\1",sources$gene_id)
    if(length(unique(nn))>1)    return(.aggSources(sources,addStatus="ambiguous"))
    w <- grep("-(3|5)p$",sources$transcript_id)
    if(length(w)>0){
        if(length(unique(sources$gene_id[w]))==1)  return(.aggSources(sources[w,,drop=F],addStatus="unique"))
        sources <- sources[w,,drop=F]
    }
    sources$gene_id <- nn
    return(.aggSources(sources,addStatus="unique"))
}


.aggSources <- function(sources, maxN=5, bypassType=NULL, addStatus=NULL){
    if(length(unique(sources$chrPos))>1) sources$chrPos <- "ambiguous"
    sources$cigar <- .agCigar(sources$cigar)
    sources$overlap <- max(sources$overlap,na.rm=T)
    sources$percentOverlap <- max(sources$percentOverlap,na.rm=T)
    sources$posInFeature <- ifelse(length(unique(sources$posInFeature))==1,sources$posInFeature[[1]],NA_integer_)
    sources <- aggregate(sources[,-1],by=list(seq=sources$seq),maxN=maxN,FUN=function(x,maxN){
        if(is.character(x)){
            x <- unique(x)
            if(length(x)>maxN) return("ambiguous")
            return(paste(sort(unique(x)),collapse=";"))
        }
        return(x[[1]])
    })
    if(!is.null(bypassType)) sources$transcript_type <- bypassType
    if(!is.null(addStatus)) sources$status <- addStatus
    return(sources)
}


# src should be a single-row data.frame of sources
.validateFeatureType <- function(src, rules){
    if(src[["transcript_type"]] %in% names(rules$typeValidation)) return(rules$typeValidation[[src[["transcript_type"]]]](src))
    return(TRUE)
}

.nmax <- function(x){ ifelse(all(is.na(x)),NA_integer_,max(x,na.rm=T)) }

.agCigar <- function(x){
    if(length(unique(x))==1) return(x[1])
    m <- sapply(x, FUN=function(x){
        x <- .splitCigar(x)
        sum(as.numeric(x[which(x[,1]=="M"),2]))
    })
    return(x[order(m,decreasing=T)[1]])
}

.disambiguate_tRNA <- function(sources, rules=defaultAssignRules()){
    if(rules$tRF$gtRNAdb=="require" & any(sources$transcript_type=="pseudo_tRNA")){
        return(.aggSources(sources,addStatus="ambiguous"))
    }else{
        if(rules$tRF$gtRNAdb=="prioritize" & any(sources$transcript_type=="tRNA")){
            sources <- sources[which(sources$transcript_type=="tRNA"),,drop=F]
            if(nrow(sources)==1){
                sources$status <- "unique"
                return(sources)
            }
        }
    }
    sources$gene_id <- gsub("_$","",sources$gene_id)
    sources$gene_id <- gsub("(m)$","",sources$gene_id)
    n2 <- gsub("-[0-9]*$","",sources$gene_id)
    n3 <- gsub("-[0-9]*$","",n2)
    type <- tRFtype(sources, rules)
    if(grepl("3p",type)){
        pos <- max(sources$posInFeature,na.rm=T)
    }else{
        if(any(sources$posInFeature==0)){
            pos <- 0
        }else{
            pos <- median(sources$chrPos,na.rm=T)
        }
    }
    if(length(unique(n2))>1 & length(unique(n3))>1){
        n4 <- .tRNAnt(n3)
        sources <- .aggSources(sources, 
                        bypassType=type,
                        addStatus=ifelse(length(unique(n4))==1,"unique","ambiguous"))
        return(sources)
    }
    if(length(unique(n2))==1){
        w <- grep("-[0-9]*-[0-9]*$",sources$gene_id)
        if(length(unique(sources$gene_id[w]))==1){
            sources <- .aggSources(sources[w,,drop=F], bypassType=type, addStatus="unique")
            if(type=="tRNA") sources$transcript_type <- tRFtype(sources[w[1],,drop=F], rules)
        }else{
            sources <- .aggSources(sources, bypassType=type, addStatus="unique")
        }
        return(sources)
    }
    w <- grep("-[0-9]*$",n2)
    if(length(unique(n2[w]))==1){
        return(.aggSources(sources[w,,drop=F], bypassType=type, addStatus="unique"))
    }
    return(.aggSources(sources, addStatus="ambiguous"))
}

.tRNAnt <- function(gnames){
    x <- strsplit(gnames,"-",fixed=T)
    x <- lapply(x, FUN=function(x){ x[[2]] <- gsub("SeC(e)","SelCys",x[[2]],fixed=T); x })
    x2 <- sapply(x, FUN=function(x){ if(length(x)>2) x[[3]] <- gsub("R","A",x[[3]],fixed=T); paste(x,collapse="-") })
    if(length(unique(x2))==1) return(x2)
    x2 <- sapply(x, FUN=function(x){ if(length(x)>2) x[[3]] <- gsub("R","G",x[[3]],fixed=T); paste(x,collapse="-") })
    if(length(unique(x2))==1) return(x2)
    x2 <- sapply(x, FUN=function(x){ if(length(x)>2) x[[3]] <- gsub("Y","C",x[[3]],fixed=T); paste(x,collapse="-") })
    if(length(unique(x2))==1) return(x2)
    x2 <- sapply(x, FUN=function(x){ if(length(x)>2) x[[3]] <- gsub("Y","T",x[[3]],fixed=T); paste(x,collapse="-") })
    if(length(unique(x2))==1) return(x2)
    return(gnames)
}

.tRNAbasename <- function(x){
    x <- gsub("-[0-9]*$","",x)
    gsub("-[0-9]*$","",x)
}
.miRNAbasename <- function(x,family=FALSE){
    x <- gsub("p(A|C|G|T)$","p",x)
    x <- gsub("-(3|5)p$","",x)
    if(family) x <- gsub("(mir-[a-z0-9]*)-[0-9]*","\\1",x)
    return(x)
}


#' checkPiRNA
#'
#' Checks whether a sequence could in principle be a piRNA
#'
#' @param seqs A character vector of sequences
#' @param type Type of sequence signature to check for. Either "primary" (default), "seconday" or "any".
#' @param allowRevComp Logical; whether to allow a sequence to match the signature with its reverse complement.
#' @param size Integer vector indicating the possible sizes (default 26:32).
#' 
#' @return A logical vector
#'
#' @export
checkPiRNA <- function(seqs, type="primary", allowRevComp=FALSE, size=26:32){
    seqs <- as.character(seqs)
    type <- match.arg(type, c("primary","secondary","any"))
    sl <- sapply(seqs, FUN=nchar)
    if(allowRevComp){
        library(Biostrings, quietly=TRUE)
        revcomp <- as.character(reverseComplement(DNAStringSet(seqs)))
        if(type=="primary" | type=="any"){
            pp <- apply(cbind(seqs,revcomp), 1, FUN=function(x){ any(sapply(strsplit(x,""),FUN=function(x){x[[1]]})=="T") })
        }
        if(type=="secondary" | type=="any"){
            sp <- apply(cbind(seqs,revcomp), 1, FUN=function(x){ any(sapply(strsplit(x,""),FUN=function(x){x[[10]]})=="A") })
        }
    }else{
        if(type=="primary" | type=="any"){
            pp <- sapply(strsplit(seqs,""),FUN=function(x){x[[1]]=="T"})
        }
        if(type=="secondary" | type=="any"){
            sp <- sapply(strsplit(seqs,""),FUN=function(x){x[[10]]=="A"})
        }
    }
    return(switch(type,
        any=( sl %in% size & (pp | sp) ),
        primary=( sl %in% size & pp ),
        secondary=( sl %in% size & sp ) ))
}


#' defaultAssignRules
#' 
#' @return A list
#'
#' @export
defaultAssignRules <- function(){
    return(list(    overlapBy=0.5,
                    sameStrand="prioritize", # either 'require', 'prioritize', or 'any'
                    prioritizeKnown=FALSE,
                    typeValidation=list(
                        primary_piRNA=function(x){ checkPiRNA(x$seq, type="primary") },
                        secondary_piRNA=function(x){ checkPiRNA(x$seq, type="secondary") },
                        miRNA=function(x){ check_miRNA(x, length=20:23, minOverlap=12) }
                    ),
                    tRF=list(gtRNAdb="prioritize", # either 'require', 'prioritize', or 'any'
                            types=list(
                            "tRNA_5p_half"=function(x){ x$posInFeature %in% -1:1 & x$length %in% 30:34 },
                            #"tRNA_3p_half"=function(x){ x$posInFeature>=(x$feature_length-x$length-4) & x$length >= 38 & x$length <= 50 },
                            "tRNA_3p_half"=function(x){ x$posInFeature>=35 & x$length >= 38 & x$length <= 50 & grepl("CCA$",x$seq)},
                            "tRNA_5p_fragment"=function(x){ x$posInFeature <5 & x$length < 30 },
                            #"tRNA_3p_fragment"=function(x){ x$posInFeature>30 & (x$feature_length-x$posInFeature-x$length)<10  },
                            "tRNA_3p_fragment"=function(x){ x$posInFeature>35 }
                    )),
                    highPriorityTypes=c("miRNA","tRNA","tRNAp","Mt_tRNA","snRNA","snoRNA","antisense","piRNA_precursor"),
                    lowPriorityTypes=c("precursor")
        ))
}


#' getPotentialSources
#'
#' Returns the potential sources (before assignment) of a sequence of interest.
#'
#' @param o An object of class shortRNAexp.
#' @param seq A character vector including the sequence(s) of interest.
#' 
#' @return A data.frame.
#'
#' @export
getPotentialSources <- function(o, seq){
    if(!is(o, "shortRNAexp")) stop("`o` should be an object of class `shortRNAexp`.")
    if(is.null(o@allsrcs) | nrow(o@allsrcs)==0) stop("The shortRNAexp object was created with keepAllSrcs=FALSE, so that only the final assigned sources were saved.")
    o@allsrcs[which(o@allsrcs$seq %in% seq),,drop=F]
}

#' tRFtype
#'
#' Returns the tRNA fragment type from a sources data.frame
#'
#' @param o An object of class shortRNAexp.
#' @param seq A character vector including the sequence(s) of interest.
#' 
#' @return A data.frame.
#'
#' @export
tRFtype <- function(srcs, rules=defaultAssignRules()){
    if(rules$tRF$gtRNAdb=="require"){
        if(length(unique(srcs$transcript_type))>1)  return("tRNA;pseudo_tRNA")
    }
    srcs2 <- srcs[grep("^tRNA",srcs$location),,drop=F]
    if(nrow(srcs2)>0){
        srcs <- srcs2
    }else{
        if(rules$tRF$gtRNAdb=="require")    return("tRNA")
    }
    a <- rep(NA_character_,nrow(srcs))
    for(i in 1:nrow(srcs)){
        for(f in names(rules$tRF$types)){
            if(rules$tRF$types[[f]](srcs[i,,drop=F])){
                a[i] <- f
            }
        }
    }
    a <- unique(a[which(!is.na(a))])
    if(length(a)==0) return("tRNA_internal_fragment")
    if(length(a)==1) return(a)
    if(all(grepl("5p",a))) return("tRNA_5p_fragment")
    if(all(grepl("3p",a))) return("tRNA_3p_fragment")  
    return("tRNA_fragment")
}
