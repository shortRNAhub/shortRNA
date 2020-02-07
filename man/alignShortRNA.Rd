% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functionsDKT.R
\name{alignShortRNA}
\alias{alignShortRNA}
\title{Alignment using Rsubread}
\usage{
alignShortRNA(
  fastq,
  fastq2 = NULL,
  index,
  nBestLocations = 100,
  nThreads = parallel::detectCores() - 1,
  GTF = NULL,
  ...
)
}
\arguments{
\item{fastq}{character vector of file names.}

\item{fastq2}{in case the sequencing is paired-end.}

\item{index}{path to the Rsubread index}

\item{nThreads}{numeric value of how many cores to be used for alignment.}

\item{GTF}{path to the GTF file.}

\item{...}{other parameters specific to `Rsubread`.}

\item{numeric}{value specifying the maximal number of equally-best mapping locations that will be reported for a multi-mapping read.}
}
\value{
Stores a `BAM` file.
}
\description{
Alignment using Rsubread
}
\seealso{
Rsubread
}
\author{
Deepak Tanwar (tanward@ethz.ch)
}