#' Dealing with soft-clipping
#'
#' @param GA `GAlignment file`
#'
#' @return `GRanges`
#' @export
#'
#' @examples
repSoftClip <- function(GA) {
  suppressPackageStartupMessages({
    library(Rsamtools)
    library(rtracklayer)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(GenomicAlignments)
    library(plyr)
    library(Biostrings)
  })
  res <- sapply(seq_along(GA), function(i) {
    # Convert GAlignment file to GRanges
    gr <- GRanges(GA[i])
    cigar <- gr$cigar
    st <- strand(gr)
    st <- as.character(st@values)
    seq <- as.character(gr$seq)
    
    if (st == "-") seq <- as.character(reverseComplement(DNAStringSet(seq)))
    
    seqCigar <- as.character(cigarToRleList(cigar)[[1]])
    nr <- seqCigar
    seqS <- strsplit(seq, "")[[1]]
    
    if (st == "-") seqS <- rev(seqS)
    
    if (length(seqCigar) != length(seqS)) {
      stop("Please provide corect sequence and matching cigar.")
    }
    
    wh <- which(seqCigar == "M")
    seqCigar[wh] <- seqS[wh]
    
    # if (st == "-") seqCigar <- rev(seqCigar)
    
    gr$seqC <- paste(seqCigar, collapse = "")
    gr$seqN <- gsub(pattern = "S", replacement = "", x = gr$seqC)
    
    nr[nr == "M"] <- start(gr):end(gr)
    
    wh1 <- which(nr != "S")
    wh1 <- wh[c(1, length(wh1))]
    se <- (start(gr) - (wh[1] - 1)):(end(gr) + (length(nr) - wh1[2]))
    start(gr) <- se[1]
    end(gr) <- se[length(se)]
    return(gr)
  })
  
  res <- GRanges(plyr::ldply(lapply(res, data.frame)))
  res
}


#' Make tracks with reads
#'
#' @param gr `GRanges` object with the genomic coordinates to be plotted
#' @param counts a `dataframe` of samples by features counts
#' @param features `GRanges` of features
#' @param bamFile path to a BAM file
#' @param plotCoverage Logical. Wether to plot coverage or not.
#'
#' @return plot with tracks
#' @export
#'
#' @examples
makeTracks <- function(gr, 
                       counts, 
                       features, 
                       bamFile = NULL, 
                       param = ScanBamParam(what = c("cigar", "pos", "seq")),
                       plotCoverage = TRUE) {
  
  suppressPackageStartupMessages({
    library(Rsamtools)
    library(ggplot2)
    library(ggbio)
    library(rtracklayer)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(GenomicAlignments)
    library(cowplot)
    library(rcartocolor)
    library(Biostrings)
  })
  
  ga <- readGAlignments(bamFile, param = param)
  names(ga) <- ga@elementMetadata$seq
  seqlevelsStyle(ga) <- "Ensembl"
  
  if (!is.null(bamFile)) bamFile <- BamFile(bamfile)
  gr <- suppressWarnings(
    subsetByOverlaps(features, reduce(gr), minoverlap = 1)
  )
  gar <- suppressWarnings(
    subsetByOverlaps(ga, gr)
  )
  
  mt <- data.frame(strand = gar@strand, gar@elementMetadata)
  mt$Counts <- NA
  
  for (i in 1:nrow(mt)) {
    if (mt$strand[i] == "+") {
      mt$Counts[i] <- mean(as.numeric(counts[mt$seq[i], ]))
    } else if (mt$strand[i] == "-") {
      mt$Counts[i] <- mean(
        as.numeric(
          counts[as.character(reverseComplement(DNAStringSet(mt$seq[i]))), ]
        )
      )
    }
  }
  
  
  gar@elementMetadata <- DataFrame(mt[, -1])
  
  gr <- repSoftClip(gar)
  
  qr <- suppressWarnings(
    subsetByOverlaps(features, reduce(gr))
  )
  qr$model <- "cds"
  names(qr) <- qr$transcript_id
  qrl <- as(qr, "GRangesList")
  
  # Make row for every letter
  df <- mapply(function(s, e, seq, i) {
    data.frame(
      pos = seq(s, e, by = 1L),
      seq = seq,
      id = i
    )
  },
  s = start(gr), e = start(gr) + nchar(gr$seq) - 1, seq = strsplit(gr$seqC, ""), i = seq_along(gr),
  SIMPLIFY = FALSE
  )
  
  # Combine different granges
  df <- do.call(rbind, df)
  
  # Determine y-position in plot
  offset <- disjointBins(gr)
  df$offset <- offset[df$id]
  
  df$seq <- factor(df$seq, levels = c("A", "T", "G", "C", "S"))
  
  # Plot
  cl <- c(RColorBrewer::brewer.pal(n = 9, name = "Set1")[c(1:4, 9)])
  names(cl) <- c("A", "T", "G", "C", "S")
  
  p1 <- autoplot(GRanges(gar),
                 geom = "arrowrect", which = qr,
                 aes(fill = Counts)
  ) +
    theme_linedraw()
  # p2 <- autoplot(gar, geom = "line", stat = "coverage", size = 1.5) +
  #   scale_color_manual(name = "Coverage", values = "black") + theme_linedraw()
  # p2
  
  p2 <- ggplot(df, aes(pos, offset, fill = seq)) +
    geom_tile(height = 0.8) +
    scale_fill_manual(values = cl) +
    theme_linedraw()
  
  if(plotCoverage){
    wh <- reduce(gr)
    seqlevelsStyle(wh) <- "UCSC"
    p3 <- suppressMessages(
      autoplot(bamFile, which = wh, size = 1, color = "black") + theme_linedraw() +
        scale_color_identity(guide = "legend")
    )  
  }
  
  
  # p4 <- ggplot() +
  #   geom_alignment(qrl, type = "model", aes(fill = transcript_type),
  #                  range.geom = "arrowrect", gap.geom = "chevron",
  #                  names.expr = "genenames:transcript_id", columns = c("genenames, transcript_id")) +
  #   scale_fill_carto_d(palette = "Safe") +
  #   theme_linedraw()
  #
  # p5 <- ggplot() +
  #   geom_alignment(qr, type = "model", aes(fill = transcript_type),
  #                  range.geom = "arrowrect", group.selfish = FALSE,
  #                  gap.geom = "chevron",
  #                  names.expr = "transcript_id") +
  #   scale_fill_carto_d(palette = "Safe") +
  #   theme_linedraw()
  
  p6 <- suppressMessages(
    autoplot(qrl, aes(type = model, fill = transcript_type, vjust = 0.1),
             range.geom = "arrowrect", label.color = "black",
             label.size = 3, hjust = 0.5, label.pos = "bottom"
    ) +
      scale_fill_carto_d(palette = "Safe") +
      theme_linedraw()
  )
  
  if(plotCoverage){
    tr <- tracks(p1, p2, p3, p6, heights = c(2, 2, 2, 2)) + ylab("")  
  } else{
    tr <- tracks(p1, p2, p6, heights = c(2, 2, 2)) + ylab("")  
  }
  
  return(tr)
}
