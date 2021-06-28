tryAdapters <- function(fq1, fq2 = NULL, adapters = NULL,
                        addReverseComplement = TRUE,
                        maxReads = 5e+05, ...) {
  library(ShortRead)
  library(stringi)
  library(Biostrings)
  library(seqTools)
  library(Rbowtie2)

  o <- identify_adapters(
    file1 = fq1, file2 = fq2,
    basename = file.path(),
    overwrite = TRUE, ...
  )


  fq <- readFastq(fq1, widthIds = F)
  if (!is.null(maxReads)) {
    fq <- fq[1:min(length(fq), maxReads), ]
  }
  if (is.null(adapters)) {
    data("adapters")
    ad <- adapters
  } else {
    if (!is(adapters, "DNAStringSet")) {
      ad <- readDNAStringSet(adapters)
    }
  }
  if (addReverseComplement) {
    ad2 <- reverseComplement(DNAStringSet(ad))
    names(ad2) <- paste(names(ad2), "(revComp)")
    ad3 <- complement(DNAStringSet(ad))
    names(ad3) <- paste(names(ad3), "(complement)")
    ad <- c(ad, ad2, ad3)
    ad <- sort(ad[which(!duplicated(as.character(ad)))])
  }
  pos <- sapply(ad, fq = fq, FUN = function(x, fq) {
    w1 <- stri_locate(fq@sread, fixed = x)[, 1]
    w2 <- width(trimLRPatterns(
      Rpattern = as.character(x),
      subject = fq, ranges = T
    ))
    apply(cbind(w1, w2), 1, na.rm = T, FUN = min)
  })
  fql <- length(fq)
  readLength <- max(width(fq)) - 1
  rm(fq)

  res <- t(apply(pos, 2, rl = readLength, FUN = function(x, rl) {
    sapply(1:rl, w = x, FUN = function(x, w) {
      sum(!(w >= x))
    })
  }))
  res <- res / fql
  colnames(res) <- 1:ncol(res)
  pos <- pos < readLength - 3
  aco <- matrix(0, ncol = length(ad), nrow = length(ad))
  diag(aco) <- 1
  colnames(aco) <- names(ad)
  rownames(aco) <- names(ad)
  for (i in 1:length(ad)) {
    for (j in 1:length(ad)) {
      if (j < i) {
        u <- sum(pos[, i] | pos[, j])
        if (u > 0) {
          aco[i, j] <- aco[j, i] <- sum(pos[, i] & pos[
            ,
            j
          ]) / u
        }
      }
    }
  }
  w <- which(res[, ncol(res) - 3] > 0.05)
  if (length(w) < 3) {
    order(res[, ncol(res) - 3], decreasing = T)[1:3]
  }
  m <- tryCatch(
    {
      library(msa)
      m <- as.data.frame(msa(ad[w], type = "dna")@unmasked)
    },
    error = function(e) {
      m <- as.data.frame(ad[w])
    }
  )
  names(m) <- "Top adapters"
  m[["% reads (>=3nt)"]] <- round(res[row.names(m), ncol(res) -
    2] * 100)
  o <- list(
    top.adapters = m, all.results = res, coocurrence = aco,
    adapters = ad, call = match.call()
  )
  class(o) <- c("adapterResults", "list")
  print(m)
  return(o)
}




plotAdapterResults <- function(o, showCoocurrence = FALSE,
                               minFreq = 0.05, addSeq = TRUE,
                               row_names_width = ifelse(showCoocurrence, 10, 16)) {
  suppressPackageStartupMessages(library(ComplexHeatmap))
  if (!is.null(minFreq)) {
    w <- row.names(o$all.results)[which(o$all.results[, ncol(o$all.results) -
      2] >= minFreq)]
  }
  if (length(w) == 0) {
    return(NULL)
  }
  if (showCoocurrence) {
    res <- o$coocurrence[w, rev(w)]
  }
  else {
    res <- o$all.results[w, ] * 100
  }
  if (addSeq) {
    row.names(res) <- paste(
      row.names(res),
      as.character(o$adapters[row.names(res)])
    )
  }
  if (showCoocurrence) {
    Heatmap(res,
      name = "Jaccard", na_col = "white",
      row_names_max_width = unit(10, "cm"),
      col = colorRampPalette(c("white", "blue", "darkblue", "black"))(30)
    )
  } else {
    cn <- 1:ncol(res)
    nbLabels <- 15
    leach <- ceiling(ncol(res) / 15)
    cn[-1 * (1:floor(length(cn) / leach)) * leach] <- ""
    colnames(res) <- cn
    Heatmap(res,
      cluster_columns = F, name = "% reads",
      row_names_max_width = unit(16, "cm"), column_title_side = "bottom",
      column_title = "Position in read",
      col = colorRampPalette(c("white", "blue", "darkblue", "black"))(30)
    )
  }
}
