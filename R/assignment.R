#' assignRead
#'
#' Assigns read to one or more of possible sources.
#'
#' @param sources A data.frame of possible sources for the read (i.e., a subset of the `sources` slot of a shortRNAexp object).
#' @param rules A list of assignment rules; defaults to `defaultAssignRules()`.
#'
#' @return A vector including the following variables: seq, location, cigar, pos_in_feature, feature_length, src name, src type, overlap length, src length, assignment status
#'
#' @export
assignRead <- function(sources, rules=defaultAssignRules()){
    sources$src_name <- gsub("mmu-miR-","mmu-mir-",as.character(sources$src_name),fixed=T)
    src_len <- sources$length[1]
    seq <- as.character(sources$sequence[1])
    if(!all(is.na(sources$location))) sources <- sources[which(!is.na(sources$location)),,drop=F]
    if(nrow(sources)==1){
        location <- sources$location
        src_name <- ifelse(sources$src_name=="unknown",NA_character_,sources$src_name)
        src_type <- ifelse(sources$src_type=="tRNA",tRFtype(sources, rules),sources$src_type)
        overlap <- sources$overlap
        pos_in_feature <- sources$pos_in_feature
        feature_length <- sources$feature_length
        cigar <- sources$cigar
    }else{
        nloc <- length(unique(sources$location))
        if(nloc>1 & nloc<nrow(sources)){
            res <- sources[1:nloc,c("sequence", "location", "cigar", "pos_in_feature", "feature_length", "src_name", "src_type", "overlap", "length")]
            res$status <- rep(NA_character_,nloc)
            for(i in 1:nloc){
                ro <- assignRead(sources[which(sources$location==unique(sources$location)[i]),,drop=F])
                for(j in 1:ncol(res))   res[i,j] <- ro[[j]]
            }
            sources <- .recastSources(res)
            sources <- res
        }
        if(any(is.na(sources$location)) | length(unique(sources$location))>1){
            if(all(is.na(sources$location))){
                location <- "unknown"
            }else{
                location <- "ambiguous"
            }
        }else{
            location <- sources$location[1]
        }
        
        if(all(is.na(sources$src_name)) | all(sources$src_name=="unknown")){
            src_name <- NA_character_
            src_type <- "unknown"
            overlap <- NA_integer_
            pos_in_feature <- NA_integer_
            feature_length <- NA_integer_
            cigar <- .agCigar(sources$cigar)
        }else{
            if(!any(is.na(sources$src_name)) & length(unique(sources$src_name))==1){
                src_name <- sources$src_name[1]
                src_type <- sources$src_type[1]
                if(src_type=="tRNA") src_type <- tRFtype(sources, rules)
                overlap <- sources$overlap[1]
                pos_in_feature <- sources$pos_in_feature[1]
                feature_length <- sources$feature_length[1]
                cigar <- .agCigar(sources$cigar)
            }else{
                res <- .disambiguateOverlappingFeatures(sources, rules=rules)
                overlap <- res[[2]]
                src_name <- res[[1]]
                src_type <- res[[3]]
                pos_in_feature <- res[[5]]
                feature_length <- res[[6]]
                cigar <- res[[4]]
            }
        }
    }
    if(!is.na(src_type) && src_type=="piRNA_precursor" && rules$primary_piRNA(seq)) src_type <- "primary_piRNA"
    if(!is.na(src_type) && src_type=="miRNA" && !.check_miRNA(sources,rules$miRNA)){
        src_type <- "ambiguous"
        status <- "ambiguous"
    }
    if(is.na(src_name)){
        if(location=="unknown" | is.na(location)){
            status <- "unmapped"
        }else{
            status <- "unknown"
        }
    }else{
        if(src_name=="ambiguous" | length(strsplit(as.character(src_name),";",fixed=T)[[1]])>1){
            status <- "ambiguous"
        }else{
            status <- "unique"
        }
    }
    return(list(seq, location, cigar, pos_in_feature, feature_length, src_name, src_type, overlap, src_len, status))
}

.disambiguateOverlappingFeatures <- function(sources, rules){
    sources <- sources[which(sources$overlap==max(sources$overlap) |
        (sources$src_type %in% rules$highPriorityTypes & sources$overlap==max(sources$overlap)-1) | 
        (sources$src_type=="miRNA" & grepl("-(3|5)p",sources$src_name) & sources$length %in% rules$miRNA$size & sources$overlap>=rules$miRNA$minOverlap) ),,drop=F]
    if(any(sources$src_type %in% rules$highPriorityTypes)){
        sources <- sources[which(sources$src_type %in% rules$highPriorityTypes),,drop=F]
    }
    dp <- unique(sources$src_name[which(sources$src_type != "precursor")])
    sources <- sources[which(!(sources$src_name %in% dp & sources$src_type=="precursor")),,drop=F]
    if(!all(sources$src_type %in% rules$lowPriorityTypes)) sources <- sources[which(!(sources$src_type %in% rules$lowPriorityTypes)),]
    if(length(unique(sources$src_name))==1){
        return(list(sources$src_name[1],sources$overlap[1],ifelse(sources$src_type[1]=="tRNA",tRFtype(sources, rules),sources$src_type[1]),.agCigar(sources$cigar),sources$pos_in_feature[1],sources$feature_length[1]))
    }
    if(all(sources$src_type == "miRNA")){
        res <- .disambiguate_miRNAs(sources, rules$miRNA$size)
        res[[5]] <- NA_integer_
        res[[6]] <- NA_integer_
        return(res)
    }
    if(all(sources$src_type %in% c("tRNAp","tRNA"))){
        return(.disambiguate_tRNA(sources, rules))
    }
    if(length(unique(sources$src_type))==1){
        src_type <- sources$src_type[1]
    }else{
        src_type <- "ambiguous"
    }
    nn <- sources$src_name
    w <- which(is.na(nn))
    if(length(w)>0) nn[w] <- sources$location[w]
    return(list(.getAmbiName(nn),.nmax(sources$overlap),src_type,.agCigar(sources$cigar),NA_integer_,NA_integer_))
}

.check_miRNA <- function(x, rules){
    if(nrow(x)>1) x <- x[1,,drop=F]
    if("size" %in% names(rules) && !(x[["length"]] %in% rules$size)) return(FALSE)
    if("minOverlap" %in% names(rules) && x[["overlap"]]<rules$minOverlap) return(FALSE)
    return(TRUE)
}

.getAmbiName <- function(nn, maxN=9){
    nn <- sort(unique(unlist(strsplit(as.character(nn),"; ",fixed=T))))
    if(length(nn)>maxN) return("ambiguous")
    paste(nn,collapse="; ")
}

.disambiguate_miRNAs <- function(sources, size=20:23){
    ov <- .nmax(sources$overlap)
    sources$src_name <- gsub("mmu-miR-","mmu-mir-",gsub("p(A|C|G|T)$","p",sources$src_name),fixed=T)
    if(!(sources$length[1] %in% size)) return(list("ambiguous",ov,"ambiguous",.agCigar(sources$cigar)))
    sources <- sources[order(sources$overlap,decreasing=T),,drop=F]
    nn <- gsub("-(3|5)p$","",sources$src_name)
    nn <- gsub("-([0-9a-b]+)-([0-9]*)$","-\\1",nn)
    if(length(unique(nn))>1) return(list(.getAmbiName(nn),ov,"miRNA",sources$cigar[1]))
    sp <- sources[grep("-(3|5)p$",sources$src_name),,drop=F]
    if(nrow(sp)==0) return(list(nn[1],ov,"miRNA",.agCigar(sources$cigar)))
    ov <- .nmax(sp$overlap)
    if(length(unique(sp$src_name))==1 | length(unique(sp$src_name[which(sp$overlap == ov)]))==1) return(list(sp$src_name[1],ov,"miRNA",.agCigar(sources$cigar)))
    return(list(nn[1],ov,"miRNA",.agCigar(sources$cigar)))
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
    sources$src_name <- gsub("_$","",sources$src_name)
    sources$src_name <- gsub("(m)$","",sources$src_name)
    sources <- sources[order(sources$overlap, decreasing=T),,drop=F]
    n2 <- gsub("-[0-9]*$","",sources$src_name)
    n3 <- gsub("-[0-9]*$","",n2)
    type <- tRFtype(sources, rules)
    if(grepl("3p",type)){
        pos <- max(sources$pos_in_feature,na.rm=T)
    }else{
        if(any(sources$pos_in_feature==0)){
            pos <- 0
        }else{
            pos <- min(sources$pos_in_feature,na.rm=T)
        }
    }
    if(length(unique(n2))>1 & length(unique(n3))>1){
        n4 <- .tRNAnt(n3)
        if(length(unique(n4))==1){
            return(list(n4[1], max(sources$overlap), type,  .agCigar(sources$cigar), pos, sources$feature_length[1]))
        }
        return(list(.getAmbiName(n3),NA_integer_, type, .agCigar(sources$cigar), NA_integer_,NA_integer_))
    }
    if(length(unique(n2))==1){
        w <- grep("-[0-9]*-[0-9]*$",sources$scr_name)
        if(length(unique(sources$scr_name[w]))==1){
            if(type=="tRNA") type <- tRFtype(sources[w[1],,drop=F], rules)
            return(list(sources$scr_name[w[1]], sources$overlap[w[1]], type, .agCigar(sources$cigar[w]), sources$pos_in_feature[w[1]], sources$feature_length[w[1]]))
        }
        return(list(n2[1], sources$overlap[1], type, .agCigar(sources$cigar), pos, sources$feature_length[1]))
    }
    w <- grep("-[0-9]*$",n2)
    if(length(unique(n2[w]))==1){
        return(list(sources$scr_name[w[1]], sources$overlap[w[1]], type, .agCigar(sources$cigar[w]), sources$pos_in_feature[w[1]], sources$feature_length[w[1]]))
    }
    return(list(n3[1], sources$overlap[1], type, .agCigar(sources$cigar), pos, sources$feature_length[1]))
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


#' aggregateSequences
#'
#' Aggregates individual sequences into features
#'
#' @param o An object of class shortRNAexp.
#' @param quiet Logical; indicates whether to suppresses messages (default FALSE).
#'
#' @return An updated object of class shortRNAexp.
#'
#' @export
aggregateSequences <- function(o, quiet=FALSE){
    if(!quiet) message("- sequences of unique origins")
    o <- .aggUnambiguous(o)

    if(!quiet) message("- sequences of ambiguous origins")
    o <- .aggAmbiguous(o)

    return(o)
}

.aggUnambiguous <- function(o){
    if(!is(o, "shortRNAexp")) stop("`o` should be an object of class `shortRNAexp`.")
    m <- o@seqcounts
    sources <- o@sources
    l1 <- .doag(m, sources, which(sources$status=="unique"))
    l2 <- .doag(m, sources, which(sources$src_type %in% c("tRNA",tRNAtype()) & sources$status=="unique"), .tRNAbasename)
    l3 <- .doag(m, sources, which(sources$src_type=="miRNA" & sources$status=="unique"), .miRNAbasename)
    l4 <- .doag(m, sources, which(sources$src_type=="miRNA" & sources$status=="unique"), function(x){ .miRNAbasename(x,T)})
    w <- which(sources$src_type %in% tRNAtype() & sources$status=="unique")
    l5 <- list( os=sources$src_name[w],
                as=paste(.tRNAbasename(sources$src_name[w]), gsub("^tRNA_","",sources$src_type[w]), sep="_"),
                ss=row.names(sources)[w]
            )
    l5$ag <- aggregate(m[l5$ss,,drop=F], by=list(l5$as), FUN=sum)
    l5$tt <- as.character(sources$src_type[w])

    ag <- aggregate(rbind(l1$ag,l2$ag,l3$ag,l4$ag,l5$ag)[,-1,drop=F], by=list(c(l1$ag[,1],l2$ag[,1],l3$ag[,1],l4$ag[,1],l5$ag[,1])), FUN=max)
    row.names(ag) <- ag[,1]
    ag[,1] <- NULL
    
    # remove duplicated rows if names are a subset
    tsig <- apply(ag,1,collapse=";",FUN=paste)
    toRemove <- unlist(sapply(unique(tsig[which(duplicated(tsig))]),tsig,ag,FUN=function(x,tsig,ag){
        e <- row.names(ag[which(tsig==x),])
        e <- e[order(nchar(e))]
        e[which(sapply(e[1:(length(e)-1)],e=e,FUN=function(y,e){ any(grepl(y,e,fixed=T)) }))]
    }))
    ag <- ag[which(!(row.names(ag) %in% toRemove)),]

    o@agcounts <- as.matrix(ag)
    
    ll <- list(l1,l2,l3,l4,l5)
    rm(l1,l2,l3,l4,l5)
    os <- unlist(lapply(ll, FUN=function(x){ x$os }))
    as <- unlist(lapply(ll, FUN=function(x){ x$as }))
    tt <- unlist(lapply(ll, FUN=function(x){ x$tt }))
    ss <- unlist(lapply(ll, FUN=function(x){ x$ss }))
    
    ag <- aggregate(cbind(os,ss,tt),by=list(as),FUN=function(x){ paste(sort(unique(as.character(x))), collapse=";") })
    row.names(ag) <- ag[,1]
    ag <- ag[row.names(o@agcounts),-1]
    colnames(ag) <- c("originalNames","sequences","type")
    o@agdef <- ag
    
    o
}

.aggAmbiguous <- function(o){
    u <- unique(o@sources$src_name[which(o@sources$status=="ambiguous" & o@sources$src_name != "ambiguous")])
    ag <- matrix(0,nrow=length(u),ncol=ncol(o@seqcounts))
    row.names(ag) <- u
    for(i in 1:length(u)){
        g <- u[i]
        sr <- o@sources[which(o@sources$src_name==g),,drop=F]
        ag[i,] <- colSums(o@seqcounts[row.names(sr),,drop=F])
        if(all(sr$src_type %in% tRNAtype())){
            type <- NULL
            if(length(unique(sr$src_type))==1){
                type <- sr$src_type[1]
            }else{
                if(all(grepl("5p",sr$src_type))){
                    type <- "5p_fragment"
                }else{
                    if(all(grepl("3p",sr$src_type))) type <- "3p_fragment"
                }
            }
            if(!is.null(type)){
                row.names(ag)[i] <- paste(g, gsub("tRNA_","",type,fixed=T))
            }
        }
    }
    o@agcounts_ambiguous <- ag

    return(o)
}

.doag <- function(m, sources, w, nFUN=function(x){x}){
    os <- sources$src_name[w]
    as <- nFUN(sources$src_name[w])
    ss <- row.names(sources)[w]
    ag <- aggregate(m[ss,,drop=F], by=list(as), FUN=sum)
    tt <- as.character(sources$src_type[w])
    return(list(os=os,as=as,ss=ss,ag=ag,tt=tt))
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
    return(list(    miRNA=list(size=20:23, minOverlap=12),
                    primary_piRNA=function(x){ checkPiRNA(x, type="primary") },
                    secondary_piRNA=function(x){ checkPiRNA(x, type="secondary") },
                    tRF=list(gtRNAdbOnly=TRUE, types=list(
                            "tRNA_5p_half"=function(x){ x$pos_in_feature<=1 & x$length %in% 30:34 },
                            "tRNA_3p_half"=function(x){ x$pos_in_feature>=(x$feature_length-x$length-4) & x$length >= 38 & x$length <= 50 },
                            "tRNA_5p_fragment"=function(x){ x$pos_in_feature<5 & x$length < 30 },
                            "tRNA_3p_fragment"=function(x){ x$pos_in_feature>30 & (x$feature_length-x$pos_in_feature-x$length)<10  }
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
#' @param srcsFile The path to the seqs.srcs file. This is needed only if the shortRNAexp object was created with keepAllSrcs=FALSE.
#' 
#' @return A data.frame.
#'
#' @export
getPotentialSources <- function(o, seq, srcsFile=NULL){
    if(!is(o, "shortRNAexp")) stop("`o` should be an object of class `shortRNAexp`.")
    if(is.null(o@allsrcs) | nrow(o@allsrcs)==0){
        if(is.null(srcsFile)) stop("The shortRNAexp object was created with keepAllSrcs=FALSE, so that only the final assigned sources were saved.")
        library(data.table)
        sr <- as.data.frame(fread(paste0('grep -w "^',seq,'" ',srcsFile)))
        if(nrow(sr)>0){
            colnames(sr) <- c("sequence","cigar","location", "pos_in_feature", "feature_length", "src_id", "src_name", "src_type", "overlap")
        }
        return(sr)
    }
    o@allsrcs[which(o@allsrcs$sequence %in% seq),,drop=F]
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
    if(rules$tRF$gtRNAdbOnly){
        srcs <- srcs[grep("^tRNA",srcs$location),,drop=F]
        if(!(nrow(srcs)>0)) return("tRNA")
    }else{
        warning("tRNA fragment identification is currently only possible for mature tRNAs (i.e. tRNAs from gtRNAdb).")
    }
    a <- rep(NA_character_,nrow(srcs))
    for(i in 1:nrow(srcs)){
        for(f in names(rules$tRF$types)){ 
            if(rules$tRF$types[[f]](srcs[i,])){
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
