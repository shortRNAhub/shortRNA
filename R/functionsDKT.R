# Build index for Rsubread from `fasta` file ----

.listFilesFTP <- function(url) {
  library(RCurl)
  library(stringr)
  files <- getURL(
    url = url,
    ftp.use.epsv = FALSE,
    dirlistonly = TRUE
  )
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
#' @param ... Other options for `Rsubread::buildindex()`.
#' 
#' @return **index** for aligning fastq files with Rsubread.
#' \dontrun{
#' @examples
#' indexRsubread(fastaGenome = "../test/reference.fa.gz", basename = "../test/reference")
#' }
#' @export
indexRsubread <- function(fastaGenome = NULL,
                          refGenome = "mm10", gencodeRelease = "M18",
                          outPath = getwd(), basename = "refRsubread",
                          ...) {
  library(Rsubread)
  library(R.utils)
  
  # If fasta file is provided
  if (!is.null(fastaGenome)) {
    message("...building index for Rsubread...\n")
    buildindex(basename = basename, reference = gunzip(fastaGenome), ...)
    gzip(gsub(pattern = "\\.gz$", "", fastaGenome))
  } else {
    # If not provided
    message(
      paste("As the `fasta` file is not provided via fastaGenome,", refGenome, "will be dowloaded and index will be build.")
      )
    
    # Various checks
    if (!grepl(pattern = "hg19|hg38|mm9|mm10", x = refGenome)) 
      stop("Please provide path to the `fasta` file! This function can only build index for mm9, mm10, hg19 or hg38, or from a `fasta file`!")
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
    download.file(url = ftp, destfile = fa)
    message("...file downloaded...\n\n...building index for Rsubread...\n")
    
    # Index building
    buildindex(basename = basename, reference = gunzip(fa), ...)
    gzip(gsub(pattern = "\\.gz$", "", fa))
  }
}


#' Read Bismark coverage files containing methylated and unmethylated read counts for CpG loci and create DGEList.
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rsubread
#'
#' @seealso Rsubread
#'
#' @param files character vector of file names.
#' @param nParallel numeric value of how many cores to be used for reading files. Important when reading files with data.table package.
#' @param verbose logical. If TRUE, read progress messages are send to standard output.
#'
#' @return A `DGEList`.
#'
#' @examples
#' dge.files <- ReadBismark2DGE(files = files, sample.names = names, data.table = T, nParallel = 8)
#'
#' @export

alignShortRNA <- function(fastq, index, ){
  nBestLocations = 1000
}