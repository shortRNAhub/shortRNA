# .RNAtoTree <- function(RNA) { library(plyr) n <- lapply(RNA$names, function(x) { a <- strsplit(x, '-')[[1]] if (a[1] ==
# 'mmu') a <- a[-1] df <- data.frame(matrix(ncol = length(a)), stringsAsFactors = F) for (i in 1:length(a)) { df[1, i] <-
# paste(a[1:i], collapse = '-') } if (length(grep(pattern = '\\.', x = df[, ncol(df)])) == 1) { b <- strsplit(df[,
# ncol(df)], '\\.')[[1]] df[, ncol(df)] <- b[1] df[, ncol(df) + 1] <- paste(b[1:2], collapse = '.') } colnames(df) <-
# c('R', paste0('C', 1:(ncol(df) - 1))) df$pathString <- paste(df[1, ], collapse = '/') return(df) }) n1 <- ldply(n,
# data.frame) return(data.frame(RNA, n1, stringsAsFactors = F)) } fastaToTree <- function(fasta) { library(Biostrings)
# library(data.tree) read.file <- readDNAStringSet(fasta) df <- data.frame(names = names(read.file), Sequence =
# as.character(read.file), stringsAsFactors = F) df$names <- as.character(sapply(df$names, function(x) strsplit(x, '\\
# ')[[1]][1])) df <- .RNAtoTree(RNA = df) tree.df <- as.Node(df) return(list(df = df, tree = tree.df)) }

# Build index for Rsubread from `fasta` file ----

.listFilesFTP <- function(url) {
  library(RCurl)
  library(stringr)
  files <- getURL(url = url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  files <- str_c(url, str_split(files, "\n")[[1]])
  files <- str_trim(files)
  return(files)
}


#' Build index for Rsubread
#' @author Deepak Tanwar (tanward@ethz.ch)
#' @import Rsubread
#' @seealso Rsubread
#' @param fastaGenome path to a gzipped fasta file (optional if want to build genome other than mm9, mm10, hg19 or hg38).
#' @param refGenome mm9/mm10/hg19/hg38. `fasta` file is downloaded from gencode (https://gencodegenes.org)
#' @param gencodeRelease Release of the genome to be downloaded.
#' @param outPath Output direcory where genome is downloaded and index is build.
#' @param basename Basename for the Rsubread index.
#' @param extraFasta `vector`. Additional fasta files location.
#' @param gzExtra if additional fasta files are gzip. Default: FALSE.
#' @param ... Other options for `Rsubread::buildindex()`.
#'
#' @return **index** for aligning fastq files with Rsubread.
#' \dontrun{
#' @examples
#' indexRsubread(fastaGenome = "../test/reference.fa.gz", basename = "../test/reference")
#' @export
indexRsubread <- function(fastaGenome = NULL, refGenome = "mm10", gencodeRelease = "M18", outPath = getwd(), basename = "refRsubread",
                          extraFasta = NULL, gzExtra = FALSE, ...) {
  library(Rsubread)
  library(R.utils)

  # If fasta file is provided
  if (!is.null(fastaGenome)) {
    message("...building index for Rsubread...\n")
    buildindex(basename = basename, reference = gunzip(fastaGenome), ...)
    gzip(gsub(pattern = "\\.gz$", "", fastaGenome))
  } else {
    # If not provided
    message(paste0(
      "As the `fasta` file is not provided via fastaGenome, ", refGenome, " will be dowloaded and index will be build from GENCODE version ",
      gencodeRelease, "."
    ))

    # Various checks
    if (!grepl(pattern = "hg19|hg38|mm9|mm10", x = refGenome)) {
      stop("Please provide path to the `fasta` file! This function can only build indexes without a `fasta` file for mm9, mm10, hg19 or hg38, or from a `fasta file`!")
    }
    if (grepl(pattern = "hg[0-9][0-9]", x = refGenome) & !grepl(pattern = "[0-9][0-9]", x = gencodeRelease)) {
      stop("Please see: https://www.gencodegenes.org/human/releases.html")
    }
    if (grepl(pattern = "mm[0-9]|mm[0-9][0-9]", x = refGenome) & !grepl(pattern = "M[0-9]|M[0-9][0-9]", x = gencodeRelease)) {
      stop("Please see: https://www.gencodegenes.org/mouse/releases.html")
    }

    # Obtaining link to download fasta file
    url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_"
    ftp <- NULL
    if (refGenome == "mm9") {
      url <- paste0(url, "mouse/release_M1/NCBIM37.genome.fa.gz")
      ftp <- url
    } else if (refGenome == "mm10") {
      u <- paste0(url, "mouse/release_", gencodeRelease, "/")
      files <- .listFilesFTP(url = u)
      ftp <- files[grep(pattern = "GRCm38.p[0-9].genome.fa.gz|GRCm38.p[0-9][0-9].genome.fa.gz", x = files)]
    } else if (refGenome == "hg19") {
      url <- paste0(url, "human/release_19/GRCh37.p13.genome.fa.gz")
      ftp <- url
    } else if (refGenome == "hg38") {
      u <- paste0(url, "human/release_", gencodeRelease, "/")
      files <- .listFilesFTP(url = u)
      ftp <- files[grep(pattern = "GRCh38.p[0-9].genome.fa.gz|GRCh38.p[0-9][0-9].genome.fa.gz", x = files)]
    }

    # Getting name to save fasta file
    str <- strsplit(x = ftp, split = "/")[[1]]
    n <- str[length(str)]
    fa <- paste0(outPath, "/", n)
    message("...downloading file...\n")
    system(paste("mkdir -p", outPath))
    download.file(url = ftp, destfile = fa)
    message("...file downloaded...\n\n...building index for Rsubread...\n")

    genomeFasta <- NULL
    # Extra fasta
    if (!is.null(extraFasta)) {
      if (gzExtra) {
        system(paste0("cat ", fa, " ", paste(extraFasta, collapse = " "), " > ", outPath, "/genome.fa.gz"))
        genomeFasta <- paste(outPath, "genome.fa.gz", sep = "/")
      } else {
        tmp <- as.character(sapply(extraFasta, gzip))
        system(paste0("cat ", fa, " ", paste(tmp, collapse = " "), " > ", outPath, "/genome.fa.gz"))
        genomeFasta <- paste(outPath, "genome.fa.gz", sep = "/")
      }
    } else {
      genomeFasta <- fa
    }

    # Index building
    buildindex(basename = paste(outPath, basename, sep = "/"), reference = genomeFasta, ...)
    # gzip(gsub(pattern = '\\.gz$', '', fa))
  }
}

# Aligning `fastq` files with Rsubread ----


#' Alignment using Rsubread
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rsubread
#'
#' @seealso Rsubread
#'
#' @param fastq character vector of file names.
#' @param fastq2 in case the sequencing is paired-end.
#' @param nThreads numeric value of how many cores to be used for alignment.
#' @param index path to the Rsubread index
#' @param numeric value specifying the maximal number of equally-best mapping locations that will be reported for a multi-mapping read.
#' @param GTF path to the GTF file.
#' @param outDir output directory. Default: current working directory
#' @param outFile basename of the output file. Default: basename of input file
#' @param mode function to be used from Rsubread package for performing alignment. Default: align
#' @param ... other parameters specific to `Rsubread`.
#' @return Stores a `BAM` file.
#'
#' @export
alignShortRNA <- function(fastq, fastq2 = NULL, index, nBestLocations = 100, nThreads = parallel::detectCores() - 1, GTF = NULL,
                          outDir = getwd(), outFile = NULL, mode = "align", ...) {
  gtfOption <- NULL
  if (!is.null(GTF)) {
    gtfOption <- TRUE
  } else {
    gtfOption <- FALSE
  }

  out <- NULL
  if (!is.null(outFile)) {
    out <- paste0(outDir, "/", outFile, ".bam")
  } else {
    o <- gsub(pattern = ".f.*", replacement = "", x = basename(fastq))
    out <- paste0(outDir, "/", o, ".bam")
  }

  library(Rsubread)
  library(parallel)

  if (mode == "align") {
    align(index = index, readfile1 = fastq, readfile2 = fastq2, nBestLocations = nBestLocations, output_file = out, nthreads = detectCores() -
      1, sortReadsByCoordinates = TRUE, annot.ext = GTF, isGTF = gtfOption, ...)
  } else if (mode == "subjunc") {
    subjunc(index = index, readfile1 = fastq, readfile2 = fastq2, nBestLocations = nBestLocations, output_file = out, nthreads = detectCores() -
      1, sortReadsByCoordinates = TRUE, annot.ext = GTF, isGTF = gtfOption, ...)
  } else {
    stop("Please provide one option: align or subjunc!")
  }
}


#' Run `lapply` function in parallel, if possible.
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import parallel
#'
#' @seealso parallel
#'
#' @param any parameters for `lapply` or `mclapply`
#' @return A `list`.
#'
#' @export
mylapply <- function(...) {
  if (require(parallel) && .Platform$OS.type == "unix") {
    mclapply(..., mc.cores = detectCores(), mc.preschedule = F)
  } else {
    lapply(...)
  }
}
