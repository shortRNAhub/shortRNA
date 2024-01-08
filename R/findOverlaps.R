#' overlapWithTx2
#'
#' Find overlaps between read alignments (bam format) and a feature annotation.
#' Supports spliced reads/transcripts, and calculates position relative to
#' transcript.
#'
#' @param bamFile A character vector of path to BAM file.
#' @param annotation An object of class `GRangesList`, as produced by
#' `prepareAnnotation()`.
#' @param ignoreStrand Logical; whether to ignore strand when searching for
#' overlaps. By default, strand is ignored and considered later on at the read
#' assignment stage.
#' @param nbthreads A positive integer indicating the number of threads to use.
#' Defaults to `min(c(8, bpworkers()))`.
#' 
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom BiocParallel SerialParam MulticoreParam bpmapply
#' @importFrom GenomicRanges setdiff intersect
#' @importFrom Rsamtools ScanBamParam
#' @importFrom ensembldb seqlevelsStyle
#' @importFrom IRanges findOverlapPairs subsetByOverlaps
#' @importFrom BiocGenerics strand start end width
#' @importFrom S4Vectors elementNROWS
#' 
#' @return A data.table.
#' @export
overlapWithTx2 <- function(bamFile, annotation,
                           ignoreStrand = TRUE, nbthreads = NULL) {
  mcols(annotation)$tx_id <- as.factor(mcols(annotation)$tx_id)

  if (is.null(nbthreads) || nbthreads == 1) {
    bp <- SerialParam()
  } else {
    bp <- MulticoreParam(nbthreads, progressbar = TRUE)
  }

  # load bam as GAlignments
  param <- ScanBamParam(what = c("cigar", "qname", "seq"))
  bam <- readGAlignments(bamFile, param = param)
  seqlevelsStyle(bam) <- "ensembl"
  bam <- as(bam, "GRangesList")
  # bam@elementMetadata$seq <- as.character(bam@elementMetadata$qname)

  message(paste(length(bam), "alignments loaded, searching for overlaps..."))

  if (is(annotation, "GRanges")) {
    annotation <- as(annotation, "GRangesList")
  }

  # suppressWarnings({
  OV <- findOverlapPairs(bam, annotation, ignore.strand = ignoreStrand)
  # OV@second <- as(OV@second, "GRangesList")
  Rdiff <- GenomicRanges::setdiff(OV)
  OVinter <- GenomicRanges::intersect(OV)
  # })

  message(paste("Found", length(OV), "overlaps."))
  message("Calculating positions relative to transcripts...")

  strands <- as.factor(unlist(unique(strand(OV@second))))
  negative_features <- which(strands == "-")
  posInFeature <- matrix(NA_integer_, ncol = 2, nrow = length(OV))
  # flag spliced alignments
  toResolve <- elementNROWS(OV@first) > 1 | elementNROWS(OV@second) > 1

  if (any(!toResolve)) {
    # unspliced
    w <- which(!toResolve)
    re <- unlist(OV@first[w, ])
    fe <- unlist(OV@second[w, ])
    posInFeature[w, 1] <- start(re) - start(fe)
    posInFeature[w, 2] <- end(fe) - end(re)
    w <- intersect(w, negative_features)
    posInFeature[w, 1:2] <- posInFeature[w, 2:1]
    rm(fe, re)
  }

  if (any(toResolve)) {
    # spliced
    w <- which(toResolve)
    posInFeature[w, ] <- t(bpmapply(
      FUN = .posInFeature,
      read = OV@first[w, ],
      feature = OV@second[w, ],
      Rd = Rdiff[w],
      oi = OVinter[w],
      BPPARAM = bp
    ))
    w <- intersect(w, negative_features)
    posInFeature[w, 1:2] <- posInFeature[w, 2:1]
  }

  rm(Rdiff)

  message("Aggregating...")

  res <- data.frame(
    seq = as.character(OV@first@elementMetadata$seq),
    cigar = OV@first@elementMetadata$cigar,
    chr = as.factor(unlist(seqnames(OV@first))),
    read.start = min(start(OV@first)),
    read.end = max(end(OV@first)),
    read.strand = as.factor(unlist(strand(OV@first))),
    overlap = sum(width(OVinter)),
    startInFeature = posInFeature[, 1],
    distanceToFeatureEnd = posInFeature[, 2],
    transcript_id = mcols(OV@second)$tx_id,
    transcript_type = factor(mcols(OV@second)$tx_biotype),
    transcript.strand = strands,
    transcript.length = sum(width(OV@second))
  )
  rm(OV, OVinter)

  # add alignments that did not overlap anything
  nonOV <- suppressWarnings(subsetByOverlaps(bam,
    annotation,
    invert = TRUE,
    ignore.strand = ignoreStrand
  ))

  
  res2 <- data.frame(
    seq = as.character(nonOV@elementMetadata$seq),
    cigar = nonOV@elementMetadata$cigar,
    chr = as.factor(unlist(unique(seqnames(nonOV)))),
    read.start = min(start(nonOV)),
    read.end = max(end(nonOV)),
    read.strand = as.factor(unlist(strand(nonOV)))
  )
  
  if(nrow(res2) > 0){
    for (f in setdiff(names(res), names(res2))) {
      res2[[f]] <- NA
    }
    res <- rbind(res, res2[, colnames(res)])
  }
  
  
  res$seq <- as.factor(res$seq)
  res
}

.posInFeature <- function(read, feature, Rd, oi) {
  c(
    .startInFeature(read, feature, Rd, oi),
    .endInFeature(read, feature, Rd, oi)
  )
}

.startInFeature <- function(read, feature, Rd, oi) {
  feature <- sort(feature)
  rs <- min(start(read))
  if (sum(width(Rd)) > 0 && rs != min(start(oi))) {
    # read start is not in transcript
    os <- rs - min(start(feature))
    if (os < 0) {
      # read starts before transcript
      return(as.integer(os))
    } else {
      # read starts in an intron
      return(NA_integer_)
    }
  }
  # read start is within transcript
  w <- which(end(feature) >= rs)[1]
  return(as.integer(sum(width(feature)[which(end(feature) < rs)]) +
    rs - start(feature)[w]))
}

.endInFeature <- function(read, feature, Rd, oi) {
  feature <- sort(feature)
  rs <- max(end(read))
  if (sum(width(Rd)) > 0 && rs != max(end(oi))) {
    # read end is not in transcript
    os <- rs - max(end(feature))
    if (os > 0) {
      # read ends after transcript
      return(as.integer(-os))
    } else {
      # read ends in an intron
      return(NA_integer_)
    }
  }
  # read end is within transcript
  w <- rev(which(start(feature) <= rs))[1]
  return(as.integer(-1 * (sum(width(feature)[which(start(feature) > rs)]) +
    end(feature)[w] - rs)))
}
