% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc.R
\name{readDuplicationPlot}
\alias{readDuplicationPlot}
\title{Reads duplication plot}
\usage{
readDuplicationPlot(
  json = NULL,
  df = NULL,
  dupRate = "dupRate",
  GCpercent = NULL,
  dupLevel = "dupLevel"
)
}
\arguments{
\item{json}{A JSON file with QC information from \code{\link{qcFastq}}}

\item{df}{A data frame}

\item{dupRate}{Column with duplication rate in df}

\item{GCpercent}{Column with GC percentage in df}

\item{dupLevel}{Column with duplication levels in df}
}
\value{
An interactive plot
}
\description{
Reads duplication plot
}
\examples{
# Input
df <- data.frame(
  dupRate = rev(sort(c(70, 10, 1:5, rep(0, 13)))),
  dupLevel = 1:20
)
df$GCpercent <- c(runif(n = 5, min = 51, max = 55), rep(0, 15))

# Analysis
## Plot with a data frame
p <- readDuplicationPlot(
  df = df,
  dupRateCol = "dupRate",
  GCpercentCol = NULL,
  dupLevelCol = "dupLevel"
)

## Plot with a data frame with GC
p_GC <- readDuplicationPlot(
  df = df,
  dupRateCol = "dupRate",
  GCpercentCol = "GCpercent",
  dupLevelCol = "dupLevel"
)

## Plot with a json file
fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
outDir <- tempdir()
qc_res <- qc_fastq(file = fq1, outDir = outDir)
p_js <- readDuplicationPlot(json = qc_res$json)

# Output
p
p_GC
p_js
}
\author{
Deepak Tanwar (tanward@ethz.ch)
}
