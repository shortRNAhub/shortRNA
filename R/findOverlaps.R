#' overlapWithTx2
#'
#' Find overlaps between read alignments (bam format) and a feature annotation.
#' Supports spliced reads/transcripts, and calculates position relative to transcript.
#'
#' @param bamFile A character vector of path to BAM file.
#' @param annotation An object of class `GRangesList`, as produced by `prepareAnnotation()`.
#' @param ignoreStrand Logical; whether to ignore strand when searching for overlaps (default FALSE)
#' @param nbthreads A positive integer indicating the number of threads to use. 
#' Defaults to `min(c(8, bpworkers()))`.
#' 
#' @return A data.table.
#' @export
overlapWithTx2 <- function(bamFile, annotation, ignoreStrand=FALSE, nbthreads=NULL){
  suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(GenomicFeatures)
    library(data.table)
    library(BiocParallel)
  })
  
  bp <- .getPara(nbthreads)
  
  # convert bam to GRanges
  param <- ScanBamParam(what=c("cigar", "seq"))
  bam <- readGAlignments(bamFile, param = param)
  bam <- as(bam, "GRangesList")
  bam@elementMetadata$seq <- as.character(bam@elementMetadata$seq)
  
  message(paste(length(bam), "alignments loaded, searching for overlaps..."))
  
  suppressWarnings({
    OV <- findOverlapPairs(bam, annotation, ignore.strand=ignoreStrand)
    Rdiff <- GenomicRanges::setdiff(OV)
    OVinter <- GenomicRanges::intersect(OV)
  })
  
  message(paste("Found", length(OV), "overlaps."))
  message("Calculating positions relative to transcripts...")
  
  strands <- .ext_rlelist(OV@second, strand)
  negative_features <- which(strands=="-")
  posInFeature <- matrix( NA_integer_, ncol=2, nrow=length(OV))
  toResolve <- elementNROWS(OV@first)>1 | elementNROWS(OV@second)>1
  
  w <- which(!toResolve)
  re <- unlist(OV@first[w,])
  fe <- unlist(OV@second[w,])
  posInFeature[w,1] <- start(re)-start(fe)
  posInFeature[w,2] <- end(fe)-end(re)
  w <- intersect(w, negative_features)
  posInFeature[w,1:2] <- posInFeature[w,2:1]
  
  rm(ov1,ov2)
  
  w <- which(toResolve)
  posInFeature[w,] <- t(bpmapply( FUN=.posInFeature, 
                                  read=OV@first[w,],
                                  feature=OV@second[w,], 
                                  Rd=Rdiff[w],
                                  oi=OVinter[w],
                                  BPPARAM=bp ))
  w <- intersect(w,negative_features)
  posInFeature[w,1:2] <- posInFeature[w,2:1]
  
  rm(Rdiff)
  
  message("Aggregating...")
  
  res <- data.table( 
      seq=as.character(OV@first@elementMetadata$seq),
      cigar=OV@first@elementMetadata$cigar,
      chr=.ext_rlelist(OV@first, seqnames),
      read.start=sapply(start(OV@first), FUN=min),
      read.end=sapply(end(OV@first), FUN=max),
      read.strand=.ext_rlelist(OV@first, strand),
      overlap=sapply(width(OVinter), FUN=sum),
      startInFeature=posInFeature[,1],
      distanceToFeatureEnd=posInFeature[,2],
      transcript_id=mcols(OV@second)$transcript_id,
      transcript_type=factor(mcols(OV@second)$transcript_type),
      transcript.strand=strands,
      gene_id=mcols(OV@second)$gene_id,
      transcript.length=sapply(width(OV@second), FUN=sum)
     )
  rm(OV, OVinter)
  
  nonOV <- suppressWarnings( subsetByOverlaps( bam, 
                                               annotation, 
                                               invert=TRUE, 
                                               ignore.strand=ignoreStrand ))
  
  res2 <- data.table( 
      seq=as.character(nonOV@elementMetadata$seq),
      cigar=nonOV@elementMetadata$cigar,
      chr=.ext_rlelist(nonOV, seqnames),
      read.start=sapply(start(nonOV), FUN=min),
      read.end=sapply(end(nonOV), FUN=max),
      read.strand=.ext_rlelist(nonOV, strand)
    )

  res <- rbind(res,res2,fill=T)
  
  res
}

.posInFeature <- function(read, feature, Rd, oi){
  c( .startInFeature(read, feature, Rd, oi),
     .endInFeature(read, feature, Rd, oi) )
}

.startInFeature <- function(read, feature, Rd, oi){
  feature <- sort(feature)
  rs <- min(start(read))
  if( sum(width(Rd))>0 && rs!=min(start(oi)) ){
    # read start is not in transcript
    os <- rs-min(start(feature))
    if(os<0){
      # read starts before transcript
      return(as.integer(os))
    }else{
      # read starts in an intron
      return(NA_integer_)
    }
  }
  # read start is within transcript
  w <- which(end(feature)>=rs)[1]
  return( as.integer( sum(width(feature)[which(end(feature)<rs)]) + 
                       rs - start(feature)[w] ) )
}

.endInFeature <- function(read, feature, Rd, oi){
  feature <- sort(feature)
  rs <- max(end(read))
  if( sum(width(Rd))>0 && rs!=max(end(oi)) ){
    # read end is not in transcript
    os <- rs-max(end(feature))
    if(os>0){
      # read ends after transcript
      return(as.integer(-os))
    }else{
      # read ends in an intron
      return(NA_integer_)
    }
  }
  # read end is within transcript
  w <- rev(which(start(feature)<=rs))[1]
  return( as.integer( -1*(sum(width(feature)[which(start(feature)>rs)]) + 
                        end(feature)[w] - rs ) ) )
}

.ext_rlelist <- function(x, acc_fun){
  as.factor(unlist(acc_fun(x))[!duplicated(rep(1:length(x), elementNROWS(x)))])
}