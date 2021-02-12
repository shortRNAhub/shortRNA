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

