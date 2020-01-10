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
  l <- mylapply(files, minl = minLength, maxl = maxLength, FUN = function(x, minl, maxl) {
    x <- readDNAStringSet(x, use.names = F, format = "fastq")
    x <- as.character(x)
    nc <- sapply(x, nchar)
    x <- x[which(nc >= minl & nc <= maxl)]
    x <- table(x)
    data.frame(seq = names(x), count = as.numeric(x), stringsAsFactors = F)
  })
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


#' shortRNAexp_align
#'
#' A wrapper for the 2-steps alignment method.
#'
#' @param fasta Path to a fasta file containing the unique reads.
#' @param outputfile Path to the output bam file.
#' @param bowtie1index Path to the base of the bowtie1 index.
#' @param starindex Path to the STAR index folder.
#' @param bowtie1 Executable for bowtie1 (default 'bowtie').
#' @param star Executable for STAR (default 'STAR').
#' @param samtools Executable for samtools (default 'samtools').
#' @param m Maximum alignments for a read to be considered (default 1000).
#' @param nthreads Number of threads for alignment (default 4).
#'
#' @export
shortRNAexp_align <- function(fasta, outputfile, bowtie1index, starindex, bowtie1 = "bowtie", star = "STAR", samtools = "samtools", m = 1000, nthreads = 4) {
  of2 <- gsub(".bam", "", outputfile, fixed = T)
  cmd <- paste0(bowtie1, " -p ", nthreads, " -v 0 -S -a --best --strata -m ", m, " -f --un ", of2, ".unmapped.fasta ", bowtie1index, " ", fasta, " | ", samtools, " view -bh > ", of2, ".unsorted.bam && 
    ", samtools, " sort -@ ", nthreads, " -m 2G ", of2, ".unsorted.bam > ", of2, ".perfectMatch.bam &&  
    rm ", of2, ".unsorted.bam")
  print(cmd)
  system(cmd)
  cmd <- paste0(star, " --genomeDir ", starindex, " --runThreadN ", nthreads, " --readFilesIn ", of2, ".unmapped.fasta --alignIntronMax 1 --outFilterMultimapNmax ", m, " --outSAMattributes NH HI NM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ", of2, ".imperfect --outSAMprimaryFlag AllBestScore && 
    ", samtools, " merge -f ", outputfile, " ", of2, ".perfectMatch.bam ", of2, ".imperfect*bam && ", samtools, " index ", outputfile)
  print(cmd)
  system(cmd)
}

#' shortRNAexp_parseBam
#'
#' Extracts reads' possible sources from an alignment file.
#'
#' @param bam Path to the sorted bam file.
#' @param elements Path to the genomic elements bed file.
#' @param outputfile Path for the output file.
#' @param shell Shell executable, default 'bash'.
#' @param samtools Executable for samtools (default 'samtools').
#' @param bedtools Executable for bedtools (default 'bedtools').
#'
#' @export
shortRNAexp_parseBam <- function(bam, elements, outputfile = NULL, shell = "bash", samtools = "samtools", bedtools = "bedtools") {
  if (is.null(outputfile)) {
    f <- tempfile("srcs")
  } else {
    f <- outputfile
  }
  p <- paste0(path.package("shortRNA"), "/parseBam.sh")
  cmd <- paste0(shell, p, bam, elements, f, samtools, bedtools)
  system(cmd)
  if (is.null(outputfile)) {
    ret <- read.delim(f, header = F, stringsAsFactors = F)
    unlink(f)
    return(ret)
  }
}