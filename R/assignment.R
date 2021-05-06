#' assignReads
#'
#' @param sources A data.frame of overlaps, as produced by
#' \code{\link{overlapWithTx}}
#' @param rules Assignment rules (see \code{\link{defaultAssignmentRules}})
#'
#' @return A \code{\link[S4Vectors]{DataFrame-class}}
#' @export
assignReads <- function(sources, rules=defaultAssignRules()){
  sources <- as.data.frame(sources)
  sources$length <- nchar(as.character(sources$seq))
  sources$seq <- as.factor(sources$seq)
  sources$transcript_type <- as.factor(sources$transcript_type)
  sources$transcript_id <- as.factor(sources$transcript_id)
  levels(sources$transcript.strand) <-
    unique(c(levels(sources$transcript.strand),"*"))
  sources <- sources[order(sources$seq,sources$chr,sources$read.start),]
  # identify valid overlaps
  sources <- getOverlapValidity(sources)
  if(!is.null(rules$reclassify) && length(rules$reclassify)>0){
    sources$reclassify <- NA_character_
    for(f in names(rules$reclassify)){
      if(length(w<-which(sources$transcript_type==f))>0){
        sources$reclassify[w] <- rules$reclassify[[f]](sources[w,,drop=FALSE])
      }
    }
  }
  flags <- c("hasValid","hasOverlap","strand")
  ambiguities <- lapply(setNames(flags,flags),FUN=function(x) list())

  # set aside fragments with zero or one overlap
  dups <- unique(sources$seq[duplicated(sources$seq)])
  dups <- which(sources$seq %in% dups)
  single.src <- sources[-dups,]
  sources <- sources[dups,]
  # is there are valid overlaps, discard invalid ones
  hasValid <- unique(sources$seq[sources$valid])
  w <- which(!sources$valid & sources$seq %in% hasValid)
  if(length(w)>0){
    ambiguities$hasValid <-
      unique(sources$seq[which(sources$seq %in% unique(sources$seq[w]))])
    sources <- sources[-w,]
  }
  # if there is overlap with known features, ignore mappings with no overlap
  hasOverlap <- unique(sources$seq[!is.na(sources$transcript_id)])
  w <- which(is.na(sources$transcript_id) & sources$seq %in% hasOverlap)
  if(length(w)>0){
    ambiguities$hasOverlap <-
      unique(sources$seq[which(sources$seq %in% unique(sources$seq[w]))])
    sources <- sources[-w,]
  }
  if(!is.null(rules$sameStrand) && rules$sameStrand=="prioritize"){
    # prioritize features on the same strand as the read
    wS <- sources$strandRead==sources$strandFeature
    w <- which(wS | !(sources$seq %in% unique(sources$seq[wS])))
    if(length(w)>0){
      ambiguities$strand <-
        unique(sources$seq[which(sources$seq %in% unique(sources$seq[w]))])
      sources <- sources[-w,]
    }
  }
  # eliminate low-priority matches
  prio <- as.integer(rules$priorities[sources$transcript_type])
  prio[is.na(prio)] <- 0L
  mP <- max(splitAsList(prio,sources$seq))
  w <- which(prio < mP[sources$seq])
  if(length(w)>0){
    ambiguities$priority <-
      unique(sources$seq[which(sources$seq %in% unique(sources$seq[w]))])
    sources <- sources[-w,]
  }

  if(isTRUE(rules$prioritizeByOverlapSize)){
    mO <- max(splitAsList(sources$overlap,sources$seq))
    w <- which(sources$overlap < mO[sources$seq])
    if(length(w)>0){
      ambiguities$overlapSize <-
        unique(sources$seq[which(sources$seq %in% unique(sources$seq[w]))])
      sources <- sources[-w,]
    }
  }

  dups <- unique(sources$seq[duplicated(sources$seq)])
  dups <- which(sources$seq %in% dups)
  single.src2 <- sources[-dups,]
  single.src2$status <- 4L # resolved ambiguity
  single.src$status <- 1L + as.integer(!is.na(single.src$transcript_id)) +
    as.integer(!is.na(single.src$chr))
  sources <- sources[dups,]
  sources$status <- 5L # ambiguous

  fields <- c("cigar","chr","read.start","read.end","overlap","startInFeature",
              "transcript.length", "distanceToFeatureEnd","transcript_id",
              "transcript.strand", "transcript_type", "reclassify")
  fields <- intersect(fields, colnames(sources))

  if(nrow(sources)>0){
    nst <- lengths(unique(splitAsList(sources$strand, sources$seq)))
    sources$transcript.strand[which(nst>1L)] <- "*"
    s2 <- DataFrame(sources[!duplicated(sources$seq),])
    for(f in fields){
      s2[[f]] <- unique(splitAsList(sources[[f]], sources$seq))[s2$seq]
    }
    sources <- s2
    rm(s2)
  }
  for(f in fields){
    if(is.factor(single.src[[f]])){
      single.src[[f]] <- as(single.src[[f]], "FactorList")
      single.src2[[f]] <- as(single.src2[[f]], "FactorList")
    }else if(is.character(single.src[[f]])){
      single.src[[f]] <- as(single.src[[f]], "CharacterList")
      single.src2[[f]] <- as(single.src2[[f]], "CharacterList")
    }else{
      single.src[[f]] <- as(single.src[[f]], "IntegerList")
      single.src2[[f]] <- as(single.src2[[f]], "IntegerList")
    }
  }
  sources <- rbind(DataFrame(single.src), DataFrame(single.src2), sources)
  for(f in fields) sources[[f]] <- sources[[f]][which(!is.na(sources[[f]]))]

  stst <- c("unmapped","noFeature","unambiguous","resolvedAmbiguity","ambiguous")
  sources$status <- factor(sources$status, seq_along(stst), stst)
  a <- vapply(ambiguities, FUN.VALUE=logical(nrow(sources)), FUN=function(x){
    sources$seq %in% x
  })
  dimnames(a) <- NULL
  a <- apply(a,1,FUN=which)
  a <- relist(factor(unlist(a),seq_along(ambiguities),names(ambiguities)),
              LogicalList(a))
  sources$resolvedAmbiguities <- a

  sources$read.strand <- Rle(sources$read.strand)
  row.names(sources) <- sources$seq
  sources$length <- sources$seq <- NULL
  if(!is.null(sources$reclassify) && any(!is.na(sources$reclassify))){
    sources$reclassify <- FactorList(sources$reclassify)
  }else{
    sources$reclassify <- NULL
  }
  sources
}

#' isValidMiRNA
#'
#' @param src A data.frame of overlaps with miRNAs, as produced by
#' \code{\link{overlapWithTx}}
#' @param rules Assignment rules (see \code{\link{defaultAssignmentRules}})
#' @param length length range
#' @param minOverlap minimum overlap with the annotated mature miRNA
#' @param maxNonOverlap maximum number of nucleotides beyond the feature
#'
#' @return A logical vector of the same length as there are rows in `src`
#' @export
isValidMiRNA <- function(src, length=19:24, minOverlap=16L, maxNonOverlap=3L){
  src$length %in% length &&
    src$overlap >= minOverlap &&
    (src$length - src$overlap) <= maxNonOverlap
}

#' isPrimaryPiRNA
#'
#' @param src A data.frame of overlaps, as produced by
#' \code{\link{overlapWithTx}}
#' @param rules Assignment rules (see \code{\link{defaultAssignmentRules}})
#' @param allowRevComp Logical; whether to allow the condition to be fulfilled
#' on the reverse complement sequences
#' @param length length range
#'
#' @return A logical vector of the same length as there are rows in `src`
#' @export
isPrimaryPiRNA <- function(src, allowRevComp=FALSE, length=26:32){
  length <- as.integer(length)
  valid <- src$length >= min(length) & src$length <= max(length)
  if(length(w <- which(valid))==0) return(valid)
  seqs <- as.character(src$seq[w])
  valid[w] <- sapply(strsplit(seqs,""),FUN=function(x){ x[[1]]=="T" })
  if(allowRevComp){
    revcomp <- as.character(reverseComplement(DNAStringSet(seqs)))
    valid[w] <- valid[w] |
      sapply(strsplit(seqs,""),FUN=function(x){ x[[1]]=="T" })
  }
  valid
}

#' isSecondaryPiRNA
#'
#' @param src A data.frame of overlaps, as produced by
#' \code{\link{overlapWithTx}}
#' @param rules Assignment rules (see \code{\link{defaultAssignmentRules}})
#' @param allowRevComp Logical; whether to allow the condition to be fulfilled
#' on the reverse complement sequences
#' @param length length range
#'
#' @return A logical vector of the same length as there are rows in `src`
#' @export
isSecondaryPiRNA <- function(src, allowRevComp=FALSE, length=26:32){
  length <- as.integer(length)
  valid <- src$length >= min(length) & src$length <= max(length)
  if(length(w <- which(valid))==0) return(valid)
  seqs <- as.character(src$seq[w])
  valid[w] <- sapply(strsplit(seqs,""),FUN=function(x){ x[[10]]=="A" })
  if(allowRevComp){
    revcomp <- as.character(reverseComplement(DNAStringSet(seqs)))
    valid[w] <- valid[w] |
      sapply(strsplit(seqs,""),FUN=function(x){ x[[10]]=="A" })
  }
  valid
}

#' getOverlapValidity
#'
#' @param sources A data.frame of overlaps, as produced by
#' \code{\link{overlapWithTx}}
#' @param rules Assignment rules (see \code{\link{defaultAssignmentRules}})
#'
#' @return Data table (same format as `sources`)
#' @export
getOverlapValidity <- function(sources, rules=defaultAssignRules()){
  if(isTRUE(rules$sameStrand=="auto")){
    stra <- .strandedness(sources, rules)
    if(isTRUE(stra$revert))
      sources$read.strand <- as.factor(ifelse(sources$read.strand=="-","+","-"))
    rules$sameStrand <- stra$rule
  }
  sources$valid <- isValidOverlap(sources, rules=rules)
  types <- list()
  if(!is.null(rules$typeValidation)) types <- rules$typeValidation
  for(typ in names(types)){
    if(!is.null(types[[typ]]$fallback) &&
       length(w <- which(sources$transcript_type == typ & !sources$valid))>0){
      if(is.factor(sources$transcript_type))
        levels(sources$transcript_type) <- unique(c(sources$transcript_type,
                                                    types[[typ]]$fallback))
      sources$transcript_type[w] <- types[[typ]]$fallback
      sources$valid[w] <- isValidOverlap(sources[w,,drop=FALSE], rules)
    }
  }
  sources
}

#' isValidOverlap
#'
#' @param srcs A data.frame of overlaps, as produced by
#' \code{\link{overlapWithTx}}
#' @param rules Assignment rules (see \code{\link{defaultAssignmentRules}})
#'
#' @return A logical vector of the same length as there are rows in `srcs`
#' @export
isValidOverlap <- function(srcs, rules=defaultAssignRules){
  valid <- rep(TRUE, nrow(srcs))
  if(!is.null(rules$overlapBy))
    valid <- srcs$overlap/srcs$length >= rules$overlapBy
  if(isTRUE(rules$sameStrand=="require"))
    valid <- valid & srcs$read.strand == srcs$transcript.strand
  types <- list()
  if(!is.null(rules$typeValidation)) types <- rules$typeValidation
  for(typ in names(types)){
    w <- which(srcs$transcript_type == typ)
    if(!is.null(types[[typ]]$length)){
      types[[typ]]$length <- as.integer(types[[typ]]$length)
      valid[w] <- valid[w] & srcs[w,"length"] %in% types[[typ]]$length
    }
    if(!is.null(types[[typ]]$fun)){
      if("length" %in% names(formals(types[[typ]]$fun)) &&
         !is.null(types[[typ]]$length)){
        valid[w] <- valid[w] & types[[typ]]$fun(srcs, length=types[[typ]]$length)
      }else{
        valid[w] <- valid[w] & types[[typ]]$fun(srcs)
      }
    }
  }
  valid
}

#' defaultAssignRules
#'
#' @param rules An optional named list of rules, overwriting parts of the
#' default rules.
#'
#' @return A list
#'
#' @export
#' @examples
#' rules <- defaultAssignRules()
defaultAssignRules <- function(rules=list()){
  defrules <- list(
    overlapBy=0.5,
    prioritizeByOverlapSize=FALSE,
    sameStrand="prioritize", # 'require', 'prioritize', or 'any' or 'auto'
    # (auto will prioritize if bias > 70%, and will require if bias > 90%)
    prioritizeKnown=TRUE,
    typeValidation=list(
      primary_piRNA=list(fun=isPrimaryPiRNA, fallback="piRNA_precursor"),
      secondary_piRNA=list(fun=isPrimaryPiRNA, fallback="piRNA_precursor"),
      miRNA=list(fun=isValidMiRNA, length=19:24, fallback="miRNA_precursor")
    ),
    reclassify=list(
      "tRNA"=tRFtype,
      "pseudo_tRNA"=tRFtype
    ),
    priorities=c(
      setNames(rep(1L,9), c("miRNA","tRNA","tRNAp","Mt_tRNA", "snRNA","snoRNA",
                            "antisense","primary_piRNA","secondary_piRNA")),
      setNames(rep(-1L,9), c("precursor","long_RNA","longRNA"))
      )
  )
  for(f in names(rules)) defrules[[f]] <- rules[[f]]
  defrules
}




#' tRFtype
#'
#' Return the type of tRNA fragment
#'
#' @param srcs A data.frame of overlaps, as produced by
#' \code{\link{overlapWithTx}}
#' @param rules Named list of validity functions, in increasing order of
#' priority
#'
#' @return A factor of tRNA fragment types
#' @export
tRFtype <- function(srcs, rules=list(
  "tRNA_internal_fragment"=function(x){ TRUE },
  "tRNA_5p_fragment"=function(x){ x$startInFeature < 5L & x$length < 30L },
  "tRNA_3p_fragment"=function(x){ x$distanceToFeatureEnd < 5L & x$length < 50L },
  "tRNA_5p_half"=function(x){ x$startInFeature %in% -1:1 & x$length %in% 30:34 },
  "tRNA_3p_half"=function(x){ x$distanceToFeatureEnd %in% -1:1 & x$length >= 34L &
      x$length <= 50L & grepl("CCA$",x$seq) }
)){
  valids <- vapply(rules, FUN.VALUE=logical(nrow(srcs)), FUN=function(fn){
    fn(srcs)
  })
  factor(apply(valids, 1, FUN=function(x) max(which(x))),
         seq_len(ncol(valids)), colnames(valids))
}


.strandedness <- function(sources, rules){
  revert <- FALSE
  flog.debug("Estimating strand match ratio")
  sbias <- sum(sources$transcript.strand==sources$read.strand,na.rm=T)/
            sum(!is.na(sources$transcript.strand))
  if(sbias>0.9){
    rules$sameStrand <- 'require'
    flog.info(paste("Strand bias",round(sbias,2),"; requiring reads and ",
                    "transcripts to be on the same strand."))
  }else{
    if(sbias>0.7){
      rules$sameStrand <- 'prioritize'
      flog.info(paste("Strand bias",round(sbias,2),"; prioritizing overlaps ",
                      "on the same strand."))
    }
  }
  if(sbias<0.3){
    revert <- TRUE
    if(sbias<0.1){
      rules$sameStrand <- 'require'
      flog.info(paste("Strand bias",round(sbias,2),"; requiring reads and ",
                      "transcripts to be on different strands."))
    }else{
      rules$sameStrand <- 'prioritize'
      flog.info(paste("Strand bias",round(sbias,2),"; prioritizing overlaps ",
                      "on different strands."))
    }
  }
  return(list(rule=rules$sameStrand, revert=revert))
}
