#' fastq2SeqCountMatrix
#'
#' Collapses a set of `.fastq` (or `fastq.gz`) files to a matrix of counts per
#'  unique sequence.
#'
#' @param files A character vector containing the path to each fastq file.
#' @param minLength The minimum length for a sequence to be included (def. 15)
#' @param maxLength The maximum length for a sequence to be included (default
#'   Inf, i.e. no maximum)
#' @param discardBelowCount Discard sequences that appear less than this number
#'   of times across all samples. The default value is 2, which means that
#'   sequences that appear only once are discarded.
#' @param BPPARAM BiocParallel BPPARAM object for multithreading (default no
#'   multithreading). Note that multithreading will require more memory and is
#'   therefore not recommended for larger datasets.
#' @param format File format, either 'fasta' or 'fastq'. If omitted will be
#'   automatically detected from the first file's name.
#' @param asSparse Logical; whether to convert the matrix to sparse format (not
#'   yet entirely supported by downstream functions)
#'
#' @return a matrix, with unique sequences as rows and samples (i.e. file
#'   basenames) as column.
#'
#' @importFrom Biostrings readDNAStringSet
#' @importFrom BiocParallel bplapply MulticoreParam
#'
#' @export
fastq2SeqCountMatrix <- function(files,
                                 minLength = 15L,
                                 maxLength = 200L,
                                 discardBelowCount = 2L,
                                 BPPARAM = BiocParallel::SerialParam(progress=TRUE),
                                 format = NULL,
                                 asSparse = FALSE) {
  fe <- sapply(files, FUN = file.exists)
  if (!all(fe)) {
    stop(paste(
      "Could not find the following input file(s): \n",
      paste(files[which(!fe)], collapse = " \n ")
    ))
  }
  if(is.null(format)){
    format <- ifelse(grepl("fasta$|fasta\\.gz$|\\.fa\\.gz$|\\.fa$",
                           head(files,1)), "fasta", "fastq")
    futile.logger::flog.info(paste("Assuming", format, "format..."))
  }
  futile.logger::flog.info(paste("Reading in", length(files), format, "file(s)"))
  l <- bplapply(files,
                minl = minLength, maxl = maxLength,
                BPPARAM = BPPARAM,
                FUN = function(x, minl, maxl) {
                  x <- readDNAStringSet(x, use.names = FALSE)
                  x <- x[which(width(x)>=minl & width(x)<=maxl)]
                  x <- table(as.factor(x))
                  data.frame(seq=names(x), count=as.integer(x),
                             stringsAsFactors=FALSE)
                }
  )
  nn <- gsub("\\.fastq$|\\.fasta$|\\.fa$|\\.fq$", "",
             gsub("\\.gz$", "", basename(files)))
  names(l) <- nn
  for (i in 1:length(l)) names(l[[i]])[2] <- nn[i]
  futile.logger::flog.info("Joining counts across samples")
  df <- plyr::join_all(l, by = "seq", type = "full", match = "first")
  row.names(df) <- df$seq
  df$seq <- NULL
  df <- as.matrix(df)
  df[which(is.na(df))] <- 0L
  df <- df[which(rowSums(df, na.rm=TRUE) >= as.integer(discardBelowCount)), ]
  if(asSparse) df <- as(df, "sparseMatrix")
  df[order(row.names(df)), ]
}

#' writeUniqueFasta
#'
#' @param countMatrix A fragment count matrix, with sequences as row.names, as
#'   produced by `fastq2SeqCountMatrix()`
#' @param fasta_file The filepath to which to save the fasta file.
#' @export
writeUniqueFasta <- function(countMatrix, fasta_file="unique_fragments.fasta"){
  fa <- DNAStringSet(row.names(countMatrix))
  names(fa) <- paste0("S", 1:length(fa))
  writeXStringSet(fa, fasta_file)
}