#' Convert RNAs annotation GRanges to a Factor List
#'
#' @param gr annotation GRanges with columns: tx_id, tx_biotype, and symbol
#' @param sp name of the species
#'
#' @return
#' @export
#'
#' @examples
otherRNA2FL <- function(gr, sp = "mm10") {
  gr <- unique(gr)
  if (all(gr$tx_id == gr$symbol)) {
    fl <- FactorList(as.list(paste(sp, gr$tx_biotype, gr$symbol)))
  } else {
    fl <- FactorList(as.list(paste(sp, gr$tx_biotype, gr$symbol, gr$tx_id)))
  }
  return(fl)
}

#' Convert tRNAs annotation GRanges to a Factor List
#'
#' @param gr annotation GRanges with columns: tx_id, tx_biotype, and symbol
#' @param sp name of the species
#'
#' @return
#' @export
#'
#' @examples
tRNAtoFL <- function(gr, sp = "mm10") {
  gr$tx_biotype <- "tRNA"
  x <- gr$symbol
  y <- lapply(x, FUN = function(x) {
    x <- strsplit(gsub(".", "-", x, fixed = TRUE), "-", fixed = TRUE)[[1]]
    x <- sapply(2:length(x), FUN = function(i) paste(x[1:i], collapse = "-"))
    return(c(sp, gr$tx_biotype[1], x))
  })
  fl <- FactorList(y)
  return(fl)
}


#' Convert miRNAs annotation GRanges to a Factor List
#'
#' @param gr annotation GRanges with columns: tx_id, tx_biotype, and symbol
#' @param sp name of the species
#'
#' @return
#' @export
#'
#' @examples
miRNAtoFL <- function(gr, spMore = TRUE, sp) {
  gr$tx_biotype <- "miRNA"
  x <- paste(gr$symbol, gr$tx_id, sep = "/")

  y <- lapply(x, FUN = function(x) {
    x1 <- strsplit(x, "/")[[1]][2]
    x2 <- strsplit(x, "/")[[1]][1]
    x <- strsplit(gsub(".", "-", x1, fixed = TRUE), "-", fixed = TRUE)[[1]]
    if (length(x) > 1) {
      if (grepl("[a-z]", x[2])) {
        x1 <- strsplit(x[2], "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl = TRUE)[[1]]
        x1 <- c(x1[1], paste0(x1[1], x1[2]))
        x <- c(x[1], x1, ifelse(test = length(x) > 2, yes = x[3:length(x)], no = ""))
        x <- x[!x %in% ""]
      }
      # x1 <- sapply(1:2, FUN = function(i) paste(x[1:i], collapse = "-"))
      x1 <- paste(x[1], x[2], sep = "-")
      if (length(x) > 2) {
        x1 <- c(x1, sapply(3:length(x), FUN = function(i) paste(x[c(1, 3:i)], collapse = "-")))
      }
      return(c(sp, gr$tx_biotype[1], x2, x1))
    } else {
      return(c(sp, gr$tx_biotype[1], x2, x))
    }
  })

  fl <- FactorList(y)
  return(fl)
}


#' Convert the annotation GRanges (obtained from `prepareAnnotation` function)
#'  to a Factor List
#' @author Deepak Tanwar (tanward@ethz.ch)
#' @param anno GRanges annotation obtained from `prepareAnnotation()`
#' @param species name for the root of the tree
#' @return A `FactorList`
#' @export
featuresAnnoToFL <- function(anno, species = "mm10") {
  if (is(anno, "GRangesList")) {
    anno <- unlist(anno)
  }
  if (!is(anno, "GRanges")) {
    stop("Please check the input.")
  }

  fl <- FactorList()

  if (any(grepl("tRNA", anno$tx_biotype))) {
    tRNA <- anno[grep("tRNA", anno$tx_biotype)]
    tRNA <- tRNAtoFL(tRNA, sp = species)
    fl <- c(fl, tRNA)
  } else if (any(grepl("miRNA", anno$tx_biotype))) {
    miRNA <- anno[grep("miRNA", anno$tx_biotype)]
    miRNA <- miRNAtoFL(miRNA, spMore = TRUE, sp = species)
    fl <- c(fl, miRNA)
  } else {
    others <- anno[grep("miRNA|tRNA", anno$tx_biotype, invert = TRUE)]
    others <- otherRNA2FL(others, sp = species)
    fl <- c(fl, others)
  }

  return(fl)
}


.objToString <- function(obj) {
  deparse(substitute(obj))
}
