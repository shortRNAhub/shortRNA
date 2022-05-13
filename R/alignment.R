#' fastq2SeqCountMatrix
#'
#' Collapses a set of `.fastq` (or `fastq.gz`) files to a matrix of counts per
#'  unique sequence.
#'
#' @param files A character vector containing the path to each fastq file.
#' @param minLength The minimum length for a sequence to be included (def. 15)
#' @param maxLength The maximum length for a sequence to be included (default
#'  Inf, i.e. no maximum)
#' @param discardBelowCount Discard sequences that appear less than this number
#' of times across all samples. The default value is 2, which means that
#' sequences that appear only once are discarded.
#'
#' @return a matrix, with unique sequences as rows and samples (i.e. file
#' basenames) as column.
#'
#' @export
fastq2SeqCountMatrix <- function(files,
                                 minLength = 15,
                                 maxLength = Inf,
                                 discardBelowCount = 2) {
  library(Biostrings)
  library(plyr)
  fe <- sapply(files, FUN = file.exists)
  if (!all(fe)) {
    stop(paste(
      "Could not find the following input file(s): \n",
      paste(files[which(!fe)], collapse = " \n ")
    ))
  }
  l <- bplapply(files,
    minl = minLength, maxl = maxLength,
    BPPARAM = MulticoreParam(8),
    FUN = function(x, minl, maxl) {
      x <- readDNAStringSet(x, use.names = F, format = "fastq")
      x <- as.character(x)
      nc <- sapply(x, nchar)
      x <- x[which(nc >= minl & nc <= maxl)]
      x <- table(x)
      data.frame(seq = names(x), count = as.numeric(x), stringsAsFactors = F)
    }
  )
  nn <- gsub("\\.fastq$", "", gsub("\\.gz$", "", basename(files)))
  names(l) <- nn
  for (i in 1:length(l)) names(l[[i]])[2] <- nn[i]
  df <- join_all(l, by = "seq", type = "full", match = "first")
  row.names(df) <- df$seq
  df$seq <- NULL
  df <- as.matrix(df)
  df[which(is.na(df))] <- 0
  df <- df[which(rowSums(df, na.rm = T) >= discardBelowCount), ]
  df[order(row.names(df)), ]
}


#' List files on a FTP server
#' @import RCurl stringr
#'
#' @param url A url of FTP location
#'
#' @return A list of files
#'
#' @examples
#' # Input
#' url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/"
#'
#' # Analysis
#' files <- listFilesFTP(url)
#'
#' # Output
#' files
listFilesFTP <- function(url) {
  library(RCurl)
  library(XML)
  # library(stringr)
  files <- getHTMLLinks(getURL(
    url = url,
    ftp.use.epsv = FALSE,
    dirlistonly = TRUE
  ))

  files <- paste0(url, files)
  return(files)
}
