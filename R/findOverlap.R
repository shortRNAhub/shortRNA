#' findoverlaps.bam.featureAnnotation
#'
#' Map a BAM file to the features file
#'
#' @param bamFile A character vector of path to BAM file.
#' @param featureAnnotation A GRanges object as produced by prepareAnnotation function.
#' @param overlapBy Minimum proportion of the read inside the overlapping feature, ranging between 0 and 1 (default: 0.5).
#' 
#' @return A dataframe.
#' @export



# Function to calculate overlap -----

findoverlaps.bam.featureAnnotation <- function(bamFile, featureAnnotation, overlapBy = 0.5){
  
  suppressPackageStartupMessages(library(GenomicAlignments))
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(Rsamtools))
  
  # convert bam to GRanges
  param <- ScanBamParam(what=c("cigar", "pos", "seq"))
  b <- readGAlignments(bamFile, param = param)
  bam <- GRanges(b)
  
  # Finding overlaps
  overlap <- suppressWarnings(findOverlaps(bam, featureAnnotation))
  read <- data.frame(bam[overlap@from,], stringsAsFactors = F)
  colnames(read) <- paste0(colnames(read), "Read")
  feature <- data.frame(featureAnnotation[overlap@to], stringsAsFactors = F)
  colnames(feature) <- paste0(colnames(feature), "Feature")
  
  m <- cbind(read, feature)
  
  # Position of the featureAnnotation
  m$posInFeature <- apply(m, 1, function(x) ifelse(x[13] == "+", 
                                                   as.numeric(x['startRead']) - as.numeric(x['startFeature']), 
                                                   as.numeric(x['endFeature']) - as.numeric(x['endRead'])))
  
  m$overlap <- apply(m, 1, function(x) min(as.numeric(x['endRead']), as.numeric(x['endFeature'])) - 
                       max(as.numeric(x['startRead']), as.numeric(x['startFeature'])) + 1)
  
  m$percentOverlap <- apply(m, 1, function(x) as.numeric(x['overlap']) / 
                              (as.numeric(x['endFeature']) - as.numeric(x['startFeature']) + 1))
  
  m <- m[m$percentOverlap >= overlapBy,]
  m$chrPos <- paste(m[,1], m[,7], sep = ":") 
  
  columns <- c("seqRead", "cigarRead", "chrPos", "strandRead", "strandFeature", "posInFeature", 
               "overlap", "percentOverlap", "transcript_idFeature", "gene_idFeature",
               "transcript_typeFeature")
  r <- m[,columns]
  r <- r[order(r[,"seqRead"]),]
  return(r)
}