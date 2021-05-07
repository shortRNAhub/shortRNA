# Running QC on FASTQ file ----

#' QC of fastq files
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#'
#' @param file character vector of file names
#' @param name character vector of output file names.
#' @param outDir output directory. Default: current working directory
#' @param adapterSeq Adapter sequence for read1
#' @param adaptersFa Adapter sequence fasta file
#' @param ... other parameters specific to `Rfastp`
#' @return A list of results and location to the JSON results file
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#'
#' # Output
#' qc_res
#' @export
qc_fastq <- function(file,
                     outDir = ".",
                     name = gsub(
                       pattern = "\\.fq|\\.fq.gz|\\.fastq|\\.fastq.gz",
                       replacement = "", x = basename(file)
                     ),
                     adapterSeq = "auto",
                     adaptersFa = "",
                     minLength = 15,
                     nThread = parallel::detectCores() - 1, ...) {
  if (adapterSeq == "auto") {
    message(
      "The adapter sequence is automatically recognised.
    If the adapter sequence is already known, it is better to specify it.\n\n"
    )
  }

  library(Rfastp)
  res <- rfastp(
    read1 = file,
    outputFastq = paste(outDir, name, sep = "/"),
    thread = nThread,
    cutMeanQual = 20,
    cutFrontMeanQual = 20,
    cutTailMeanQual = 20,
    cutLowQualTail = T,
    cutLowQualFront = T,
    qualityFiltering = T,
    maxNfilter = 5,
    qualityFilterPhred = 20,
    qualityFilterPercent = 40,
    lengthFiltering = T,
    minReadLength = minLength,
    maxIndexMismatch = 1,
    adapterTrimming = TRUE,
    adapterSequenceRead1 = adapterSeq,
    adapterFasta = adaptersFa,
    ...
  )

  json <- list.files(path = outDir, pattern = name, full.names = T)
  json <- grep(pattern = "json", x = json, value = TRUE)

  ret <- list(
    result = res,
    json = json
  )
  return(ret)
}



# Various tables with some QC values ----

#' Make general information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#'
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#' df <- make_general_df(json = qc_res$json)
#'
#' # Output
#' df
#' @export
make_general_df <- function(json) {
  library(rjson)
  js <- fromJSON(file = json)

  # pkg <- as.character(packageVersion("Rfastp"))
  seq <- js$read1_before_filtering$total_cycles
  seq_type <- ifelse(test = any(grepl(pattern = "read2", x = names(js))),
    yes = "paired-end", no = "single-end"
  )
  cycle <- paste0(seq_type, " (", seq, " cycles)")

  mean_len_bf <- js$summary$before_filtering$read1_mean_length
  mean_len_af <- js$summary$after_filtering$read1_mean_length
  dup_rate <- ifelse(test = seq_type == "paired-end",
    yes = paste0(round(x = (js$duplication$rate * 100), digits = 2), "%"),
    no = paste0(
      round(x = (js$duplication$rate * 100), digits = 2),
      "% (may be overestimated since this is SE data)"
    )
  )

  df <- data.frame(
    test = c(
      # "Rfastp version",
      "Sequencing",
      "Mean length before filtering",
      "Mean length after filtering",
      "Duplication rate"
    ),
    value = c(
      # pkg,
      cycle,
      mean_len_bf,
      mean_len_af,
      dup_rate
    )
  )
  return(df)
}

# make_general_df(json = qc_res$json)


#' Make before QC information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#'
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#' df <- make_before_filt_df(json = qc_res$json)
#'
#' # Output
#' df
#' @export
make_before_filt_df <- function(json) {
  js <- rjson::fromJSON(file = json)
  tot_reads <- js$read1_before_filtering$total_reads
  tot_bases <- js$read1_before_filtering$total_bases
  q_20_bases <- js$read1_before_filtering$q20_bases
  q_20_rate <- js$summary$before_filtering$q20_rate * 100
  q_30_bases <- js$read1_before_filtering$q30_bases
  q_30_rate <- js$summary$before_filtering$q30_rate * 100
  gc <- paste0(round(x = js$summary$before_filtering$gc_content * 100, digits = 2), "%")

  df <- data.frame(
    test = c(
      "Total reads",
      "Total bases",
      "Q20 bases",
      "Q30 bases",
      "GC content"
    ),
    value = c(
      tot_reads,
      tot_bases,
      paste0(q_20_bases, " (", q_20_rate, "%)"),
      paste0(q_30_bases, " (", q_30_rate, "%)"),
      gc
    )
  )
  return(df)
}

# make_before_filt_df(json = qc_res$json)


#' Make after QC information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#' df <- make_after_filt_df(json = qc_res$json)
#'
#' # Output
#' df
#' @export
make_after_filt_df <- function(json) {
  js <- rjson::fromJSON(file = json)
  tot_reads <- js$read1_after_filtering$total_reads
  tot_bases <- js$read1_after_filtering$total_bases
  q_20_bases <- js$read1_after_filtering$q20_bases
  q_20_rate <- js$summary$after_filtering$q20_rate * 100
  q_30_bases <- js$read1_after_filtering$q30_bases
  q_30_rate <- js$summary$after_filtering$q30_rate * 100
  gc <- paste0(
    round(x = js$summary$after_filtering$gc_content * 100, digits = 2),
    "%"
  )

  df <- data.frame(
    test = c(
      "Total reads",
      "Total bases",
      "Q20 bases",
      "Q30 bases",
      "GC content"
    ),
    value = c(
      tot_reads,
      tot_bases,
      paste0(q_20_bases, " (", q_20_rate, "%)"),
      paste0(q_30_bases, " (", q_30_rate, "%)"),
      gc
    )
  )
  return(df)
}

# make_after_filt_df(json = qc_res$json)


#' Make filtering information information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#' df <- make_filt_df(json = qc_res$json)
#'
#' # Output
#' df
#' @export
make_filt_df <- function(json) {
  js <- rjson::fromJSON(file = json)

  reads <- js$summary$before_filtering$total_reads
  read_pass <- js$filtering_result$passed_filter_reads
  read_pass_pct <- round((read_pass / reads) * 100, 2)
  read_low_Q <- js$filtering_result$low_quality_reads
  read_low_Q_pct <- round((read_low_Q / reads) * 100, 2)
  read_with_N <- js$filtering_result$too_many_N_reads
  read_with_N_pct <- round((read_with_N / reads) * 100, 2)
  read_short <- js$filtering_result$too_short_reads
  read_short_pct <- round((read_short / reads) * 100, 2)
  read_long <- js$filtering_result$too_long_reads
  read_long_pct <- round((read_pass / reads) * 100, 2)

  df <- data.frame(
    test = c(
      "Reads passed filters",
      "Reads with low quality",
      "Reads with too many N",
      "Reads too short",
      "Reads too long"
    ),
    value = c(
      paste0(read_pass, " (", read_pass_pct, "%)"),
      paste0(read_low_Q, " (", read_low_Q_pct, "%)"),
      paste0(read_with_N, " (", read_with_N_pct, "%)"),
      paste0(read_short, " (", read_short_pct, "%)"),
      paste0(read_long, " (", read_long_pct, "%)")
    )
  )
  return(df)
}

# make_filt_df(json = qc_res$json)


#' Duplicated reads plot
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param df A data frame
#' @param dupRateCol Column with duplication rate in df
#' @param GCpercentCol Column with GC percentage in df
#' @param dupLevelCol Column with duplication levels in df
#'
#' @return An interactive plot
#'
#' @examples
#' # Input
#' df <- data.frame(
#'   dupRate = rev(sort(c(70, 10, 1:5, rep(0, 13)))),
#'   dupLevel = 1:20
#' )
#' df$GCpercent <- c(runif(n = 5, min = 51, max = 55), rep(0, 15))
#'
#' # Analysis
#' p <- reads_duplication_plot(
#'   df = df,
#'   dupRateCol = "dupRate",
#'   GCpercentCol = NULL,
#'   dupLevelCol = "dupLevel"
#' )
#'
#' p_GC <- reads_duplication_plot(
#'   df = df,
#'   dupRateCol = "dupRate",
#'   GCpercentCol = "GCpercent",
#'   dupLevelCol = "dupLevel"
#' )
#'
#' # Output
#' p
#' p_GC
#' @export
reads_duplication_plot <- function(df,
                                   dupRateCol,
                                   GCpercentCol = NULL,
                                   dupLevelCol) {
  library(dplyr)
  library(plotly)

  p <- plot_ly(df, x = df[, dupLevelCol]) %>%
    add_trace(
      y = df[, dupRateCol],
      type = "bar",
      name = "Read percent",
      hoverinfo = "text+x",
      hovertext = paste(
        "<b>", df[, dupRateCol], "</b>"
      )
    )

  if (!is.null(GCpercentCol)) {
    df[, GCpercentCol][df[, GCpercentCol] == 0] <- ""
    p <- p %>%
      add_lines(
        y = df[, GCpercentCol],
        name = "Mean GC ratio (%)",
        line = list(color = "red", width = 3),
        hoverinfo = "text+x",
        hovertext = paste(
          "<b>", df[, GCpercentCol], "</b>"
        )
      ) %>%
      layout(
        title = paste0(
          "Mean duplication rate: ",
          round(mean(df[, dupRateCol]), 2),
          "%"
        ),
        xaxis = list(title = "Duplication level"),
        yaxis = list(title = "Reads percent & GC ratio (%)"),
        hovermode = "x closest"
      )
  } else {
    p <- p %>%
      layout(
        title = paste0(
          "Mean duplication rate: ",
          round(mean(df[, dupRateCol]), 2),
          "%"
        ),
        xaxis = list(title = "Duplication level"),
        yaxis = list(title = "Reads percent"),
        hovermode = "x closest"
      )
  }

  return(p)
}

#' Make duplicated reads plot from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param json A JSON file
#'
#' @return An interactive plot
#'
#' #' @examples
#' # Input
#' fq1 <- system.file("extdata", "Fox3_Std_small.fq.gz", package = "Rfastp")
#' outDir <- tempdir()
#'
#' # Analysis
#' qc_res <- qc_fastq(file = fq1, outDir = outDir)
#' p <- make_dup_plot(json = qc_res$json)
#'
#' # Output
#' p
#'
#' @export
make_dup_plot <- function(json) {
  library(rjson)

  js <- fromJSON(file = json)

  bars <- js$duplication$histogram
  gc <- js$duplication$mean_gc
  xaxis <- 1:length(bars)

  data <- data.frame(
    Duplicate = round(bars / sum(bars) * 100, 2),
    GC = round(gc * 100, 2),
    x = xaxis
  )

  p <- reads_duplication_plot(
    df = data,
    dupRateCol = "Duplicate",
    GCpercentCol = "GC",
    dupLevelCol = "x"
  )
  return(
    p
    # tagList = htmltools::tagList(p)
  )
}

# a <- make_dup_plot(json = qc_res$json)



## Read quality plot
#' Make filtering information information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' \dontrun{
#' # Input
#'
#' # Analysis
#'
#' # Output
#' }
#'
#' @export
reads_quality_plot <- function(df,
                               bpPosCol,
                               A_col,
                               T_col,
                               G_col,
                               C_col,
                               meanQual = NULL) {
  library(plotly)
  library(dplyr)

  if (is.null(meanQual)) {
    df$meanQual <- rowMeans(df[, c(A_col, T_col, G_col, C_col)])
    meanQual <- "meanQual"
  }

  df$Q20 <- 20
  df$Q28 <- 8
  df$Q50 <- 18

  cols <- rcartocolor::carto_pal(n = 10, name = "Safe")[c(2, 1, 3, 4, 10)]

  plot_ly(df, x = df[, bpPosCol]) %>%
    add_trace(
      y = df[, "Q20"], name = "Not good", fillcolor = "rgba(255,0,0, 0.2)",
      type = "scatter", mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, "Q28"], name = "OK", fillcolor = "rgba(255,255,0,0.2)",
      type = "scatter", mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, "Q50"], name = "Good", fillcolor = "rgba(0,128,0,0.2)",
      type = "scatter", mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_lines(
      y = df[, A_col],
      name = "A",
      line = list(color = cols[1], width = 2),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, A_col], "</b>")
    ) %>%
    add_lines(
      y = df[, T_col],
      name = "T",
      line = list(color = cols[2], width = 2),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, T_col], "</b>")
    ) %>%
    add_lines(
      y = df[, G_col],
      name = "G",
      line = list(color = cols[3], width = 2),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, G_col], "</b>")
    ) %>%
    add_lines(
      y = df[, C_col],
      name = "C",
      line = list(color = cols[4], width = 2),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, C_col], "</b>")
    ) %>%
    add_lines(
      y = df[, meanQual],
      name = "Mean",
      line = list(color = "black", width = 3, dash = "dot"),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, meanQual], "</b>")
    ) %>%
    layout(
      title = "Per base sequence quality",
      xaxis = list(title = "Position in read (bp)"),
      yaxis = list(title = "Quality score"),
      hovermode = "x closest"
    )
}

#' Make filtering information information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' \dontrun{
#' # Input
#'
#' # Analysis
#'
#' # Output
#' }
#'
#' @export
make_read_qual_plot <- function(json, which = "before") {
  js <- rjson::fromJSON(file = json)

  qual <- NULL
  if (which == "before") {
    qual <- js$read1_before_filtering$quality_curves
  } else {
    qual <- js$read1_after_filtering$quality_curves
  }

  data <- data.frame(
    pos = 1:length(qual$mean),
    Aa = qual$A,
    Tt = qual$T,
    Gg = qual$G,
    Cc = qual$C,
    Mm = qual$mean
  )

  p <- reads_quality_plot(
    df = data,
    bpPosCol = "pos",
    A_col = "Aa",
    T_col = "Tt",
    G_col = "Gg",
    C_col = "Cc",
    meanQual = "Mm"
  )
  return(p)
}

# make_read_qual_plot(json = qc_res$json, which = "before")
# make_read_qual_plot(json = qc_res$json, which = "after")



## Base content line plot ---
#' Make filtering information information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' \dontrun{
#' # Input
#'
#' # Analysis
#'
#' # Output
#' }
#'
#' @export
base_ratio_line_plot <- function(df,
                                 bpPosCol,
                                 A_col,
                                 T_col,
                                 G_col,
                                 C_col,
                                 GC_col = NULL,
                                 N_col,
                                 isRatio = FALSE) {
  if (isRatio) {
    df <- df[, c(bpPosCol, A_col, T_col, G_col, C_col, N_col)] * 100
  }

  if (is.null(GC_col)) {
    df$GC_col <- rowMeans(df[, c(G_col, C_col)])
    GC_col <- "GC_col"
  }

  cols <- rcartocolor::carto_pal(n = 11, name = "Safe")[c(2, 1, 3, 4, 10)]

  plot_ly(df, x = df[, bpPosCol]) %>%
    add_lines(
      y = df[, A_col],
      name = "A",
      line = list(color = cols[1], width = 2),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, A_col], "</b>")
    ) %>%
    add_lines(
      y = df[, T_col],
      name = "T",
      line = list(color = cols[2], width = 2),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, T_col], "</b>")
    ) %>%
    add_lines(
      y = df[, G_col],
      name = "G",
      line = list(color = cols[3], width = 2),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, G_col], "</b>")
    ) %>%
    add_lines(
      y = df[, C_col],
      name = "C",
      line = list(color = cols[4], width = 2),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, C_col], "</b>")
    ) %>%
    add_lines(
      y = df[, GC_col],
      name = "GC",
      line = list(color = cols[5], width = 3, dash = "dot"),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, GC_col], "</b>")
    ) %>%
    add_lines(
      y = df[, N_col],
      name = "N",
      line = list(color = "black", width = 3),
      hoverinfo = "text+x",
      hovertext = paste("<b>", df[, N_col], "</b>")
    ) %>%
    layout(
      title = "Per base sequence content",
      xaxis = list(title = "Position in read (bp)/ Cycle"),
      yaxis = list(title = "Percentage (%)"),
      hovermode = "x closest"
    )
}

#' Make filtering information information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' \dontrun{
#' # Input
#'
#' # Analysis
#'
#' # Output
#' }
#'
#' @export
make_base_ratio_line_plot <- function(json, which = "before") {
  js <- rjson::fromJSON(file = json)

  qual <- NULL
  if (which == "before") {
    qual <- js$read1_before_filtering$content_curves
  } else {
    qual <- js$read1_after_filtering$content_curves
  }

  data <- data.frame(
    pos = 1:length(qual$A),
    Aa = qual$A * 100,
    Tt = qual$T * 100,
    Gg = qual$G * 100,
    Cc = qual$C * 100,
    GC = qual$GC * 100,
    Nn = qual$N * 100
  )

  p <- base_ratio_line_plot(
    df = data,
    bpPosCol = "pos",
    A_col = "Aa",
    T_col = "Tt",
    G_col = "Gg",
    C_col = "Cc",
    GC_col = "GC",
    N_col = "Nn",
    isRatio = FALSE
  )

  return(p)
}


## Base content ratio proportion plot ----
#' Make filtering information information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' \dontrun{
#' # Input
#'
#' # Analysis
#'
#' # Output
#' }
#'
#' @export
base_ratio_proportion_plot <- function(df,
                                       bpPosCol,
                                       A_col,
                                       T_col,
                                       G_col,
                                       C_col,
                                       GC_col = NULL,
                                       N_col,
                                       isRatio = FALSE) {
  if (isRatio) {
    df <- df[, c(bpPosCol, A_col, T_col, G_col, C_col, N_col)] * 100
  }

  if (is.null(GC_col)) {
    df$GC_col <- rowMeans(df[, c(G_col, C_col)])
    GC_col <- "GC_col"
  }

  cols <- rcartocolor::carto_pal(n = 11, name = "Safe")[c(2, 1, 3, 4, 10)]

  plot_ly(data, x = df[, bpPosCol]) %>%
    add_trace(
      y = df[, N_col], name = "N", fillcolor = "black", type = "scatter",
      mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, C_col], name = "C", fillcolor = cols[4], type = "scatter",
      mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, G_col], name = "G", fillcolor = cols[3], type = "scatter",
      mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, T_col], name = "T", fillcolor = cols[2], type = "scatter",
      mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_trace(
      y = df[, A_col], name = "A", fillcolor = cols[1], type = "scatter",
      mode = "none", stackgroup = "one", hoverinfo = "skip"
    ) %>%
    add_lines(
      y = df[, GC_col],
      name = "GC",
      line = list(color = cols[5], width = 3, dash = "dot"),
      hoverinfo = "text",
      hovertext = paste("<b>", df[, GC_col], "</b>")
    ) %>%
    layout(
      title = "Per base sequence content",
      xaxis = list(title = "Position in read (bp)/ Cycle"),
      yaxis = list(title = "Percentage (%)"),
      hovermode = "x closest"
    )
}


#' Make filtering information information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' \dontrun{
#' # Input
#'
#' # Analysis
#'
#' # Output
#' }
#'
#' @export
make_base_ratio_prop_plot <- function(json, which = "before") {
  js <- rjson::fromJSON(file = json)

  qual <- NULL
  if (which == "before") {
    qual <- js$read1_before_filtering$content_curves
  } else {
    qual <- js$read1_after_filtering$content_curves
  }

  data <- data.frame(
    pos = 1:length(qual$A),
    Aa = qual$A * 100,
    Tt = qual$T * 100,
    Gg = qual$G * 100,
    Cc = qual$C * 100,
    GC = qual$GC * 100,
    Nn = qual$N * 100
  )

  p <- base_ratio_proportion_plot(
    df = data,
    bpPosCol = "pos",
    A_col = "Aa",
    T_col = "Tt",
    G_col = "Gg",
    C_col = "Cc",
    GC_col = "GC",
    N_col = "Nn",
    isRatio = FALSE
  )

  return(p)
}


#' Make filtering information information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' \dontrun{
#' # Input
#'
#' # Analysis
#'
#' # Output
#' }
#'
#' @export
# `qckitfastq`
run_qc <- function(fq) {
  out <- list()

  infile <- fq
  fseq <- seqTools::fastqq(infile)

  ## Read length
  read_len <- read_length(fseq)
  # plot_read_length(read_len)

  ## Per base sequence quality
  bs <- per_base_quality(infile)
  # plot_per_base_quality(bs)

  ## Per read quality
  prq <- per_read_quality(infile)
  # plot_per_read_quality(prq)

  ## GC content
  gc_df <- GC_content(infile)
  # plot_GC_content(gc_df)

  ## Nucleotide read content
  rc <- read_content(fseq)
  # plot_read_content(rc)

  ## Kmer count
  km <- kmer_count(infile, k = 6)

  ## Overrep reads
  # overrep_reads <- overrep_reads(infile)
  # plot_overrep_reads(overrep_reads)

  ## Overrep kmer
  overkm <- overrep_kmer(infile, 7)
  # plot_overrep_kmer(overkm)

  ## Adapter content
  if (.Platform$OS.type != "windows") {
    ac_sorted <- adapter_content(infile)
    # plot_adapter_content(ac_sorted)
  } else {
    ac_sorted <- "adapter_content not available for Windows; skipping"
    print("adapter_content not available for Windows; skipping")
  }

  out <- list(
    read_length = read_len,
    per_base_quality = bs,
    per_read_quality = prq,
    GC_content = gc_df,
    read_content = rc,
    kmer_count = km,
    overrep_kmer = overkm,
    overrep_reads = overrep_reads,
    adapter_content = ac_sorted
  )
  return(out)
}







## Read length distribution
#' Make filtering information information data frame from JSON file
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param json A JSON file with QC information
#'
#' @return A data frame with FASTQ file information
#'
#' @examples
#' \dontrun{
#' # Input
#'
#' # Analysis
#'
#' # Output
#' }
#'
#' @export
read_length_hist <- function(df, groupName, readLength, nRead, col) {
  plot_ly(hoverinfo = "text", type = "bar", textposition = "auto") %>%
    add_trace(
      data = df,
      legendgroup = df[, groupName],
      x = df[, readLength],
      y = df[, nRead],
      name = df[, groupName],
      marker = list(color = col),
      hovertext = paste(
        "<b>Length:</b> ", df[, readLength],
        "<br><b>nReads:</b> ", df[, nRead]
      )
    ) %>%
    layout(
      title = "Histogram of read length distribution",
      xaxis = list(title = "Reads length"),
      yaxis = list(title = "No. of reads"),
      hovermode = "x closest"
    )
}

#' Read length distribution plot
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import Rfastp
#'
#' @seealso Rfastp
#' @param before_trim A dataframe with 2 columns: read_length and num_reads
#' @param after_trim A dataframe with 2 columns: read_length and num_reads
#'
#' @return An interactive barplot
#'
#' @examples
#' \dontrun{
#' # Input
#'
#' # Analysis
#'
#' # Output
#' }
#'
#' @export
make_read_length_plot <- function(before_trim, after_trim) {
  if (!colnames(before_trim) %in% c("read_length", "num_reads") |
    !colnames(after_trim) %in% c("read_length", "num_reads")) {
    message("Please check the columns of input")
  }
  raw <- before_trim
  trim <- after_trim

  p1 <- read_length_hist(
    df = raw,
    groupName = "Raw",
    readLength = "read_length",
    nRead = "num_reads",
    col = "red"
  )

  p2 <- read_length_hist(
    df = trim,
    groupName = "Trimmed",
    readLength = "read_length",
    nRead = "num_reads",
    col = "green"
  )

  p <- subplot(p1, p2, shareX = T, shareY = T)
  return(p)
}



## Comparison of various FastQ file reading tools "/mnt/IM/projects/software/ssc_shortRNA/04_comp_fastq_read"
