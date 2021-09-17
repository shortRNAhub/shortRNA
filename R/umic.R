#' umiCollapse
#' This function is adapted from the code available on: 
#' https://github.com/BiodataAnalysisGroup/UMIc. The function is adapted as the
#' original code was throwing some errors and was extremely slow for big `fastq`
#' files. Here, we try to make it a bit faster.
#' 
#' @param fastq path to a fastq file
#' @param UMIlocation UMI located. Default: in Read1 --> "R1"
#' @param UMIlength UMI length in base-pair. Default: 8
#' @param sequenceLength Length of sequences. Default: 51bp
#' @param countsCutoff Minimum read counts per UMI, for initial data cleaning.
#' Default: 5
#' @param UMIdistance Maximum UMI distance for UMI merging. Default: 1
#' @param sequenceDistance Maximum sequence distance for UMI correction.
#' Default: 3
#' @param outDir Output directory for results. Default: current working
#'  directory
#'
#' @return
#' @export
#'
#' @examples
#' ## Input
#' fastq <- system.file("extdata", "case3_R1.fastq.gz", package = "shortRNA")
#' 
#' ## Analysis
#' umiCollapse(fastq = fastq, UMIlength = 12, sequenceLength = 251, 
#' outDir = getwd())
#' 
umiCollapse <- function(fastq,
                        UMIlocation = "R1",
                        UMIlength = 8,
                        sequenceLength = 51,
                        countsCutoff = 5,
                        UMIdistance = 1,
                        sequenceDistance = 3,
                        outDir = getwd()) {

  # Packages required
  suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(ShortRead)
    library(Biostrings)
    library(stringdist)
    library(pryr)
    library(parallel)
  })

  message("Loaded required libraries successfully!")
  
  # Results directory and log file
  if (!dir.exists(outDir)) dir.create(outDir)
  file.create(paste0(outDir, "/extra_info.txt"))

  memory <- c()
  memory[1] <- mem_used()

  # read input file
  reads1 <- readFastq(fastq)
  
  message("FastQ file read!")

  start.time <- Sys.time()
  start.all <- proc.time()

  # File 1
  seq <- as.data.table(sread(reads1))
  ids <- as.data.table(reads1@id)

  full <- cbind(seq, ids)
  names(full) <- c("seq", "id")

  # separate UMI and read
  full$UMI <- substring(full$seq, 1, UMIlength)
  full$read <- substring(full$seq, UMIlength + 1, sequenceLength)

  full <- select(full, read, id, UMI)
  colnames(full) <- c("read", "id", "UMI")


  # keep intermediate information
  partIDs <- as.data.table(cbind(UMI = full$UMI, ID1 = full$id))
  partIDs <- partIDs[!duplicated(UMI), ]
  partIDs <- partIDs[order(partIDs$UMI, decreasing = TRUE), ]

  intermediate.table <- full[, .(count = .N), by = UMI, ]
  intermediate.table <- intermediate.table[order(intermediate.table$UMI,
    decreasing = TRUE
  ), ]
  intermediate.table$ID1 <- partIDs$ID1
  intermediate.table <- intermediate.table[
    which(intermediate.table$count >= countsCutoff),
  ]

  rm(partIDs)
  memory[2] <- mem_used()

  intermediate.table <- intermediate.table[
    order(intermediate.table$count, decreasing = TRUE),
  ]

  q <- quality(reads1)
  q1 <- split(q, c(1:10))
  q2 <- lapply(q1, as, "matrix")
  quality <- plyr::ldply(q2, as.data.table, .id = NULL)

  quality <- quality[, (UMIlength + 1):sequenceLength]
  quality$id <- full$id

  rm(ids, seq, reads1)

  # first consensus
  result_mean <- .groupingSingle(
    intermediate.table = intermediate.table,
    full = full,
    quality = quality,
    UMIlength = UMIlength
  )
  
  message("First consensus!")

  memory[3] <- mem_used()
  
  # UMI correction
  newUMIs <- .UMIcorrectionSingle(
    intermediate.table = intermediate.table,
    first_consensus = result_mean,
    sequenceDistance = sequenceDistance,
    UMIdistance = UMIdistance,
    outDir = outDir
  )

  message("UMI correction!")
  
  rm(intermediate.table)
  memory[4] <- mem_used()
  
  # final consensus
  consensus_mean <- .groupingFinalSingle(
    newUMIs = newUMIs,
    full = full,
    quality = quality,
    first_consensus = result_mean,
    UMIlength = UMIlength
  )
  
  message("Final consensus!")

  memory[5] <- mem_used()

  end.time <- Sys.time()
  end.all <- proc.time()

  total.time.secs <- difftime(end.time, start.time, units = "secs")
  total.time.mins <- difftime(end.time, start.time, units = "mins")
  time.fifference <- end.all - start.all

  line <- paste0("Total time in secs: ", total.time.secs)
  write(line, paste0(outDir, "/extra_info.txt"), append = T)

  line <- paste0("Total time in mins: ", total.time.mins)
  write(line, paste0(outDir, "/extra_info.txt"), append = T)

  line <- paste0("User time: ", time.fifference[1])
  write(line, paste0(outDir, "/extra_info.txt"), append = T)

  line <- paste0("System time: ", time.fifference[2])
  write(line, paste0(outDir, "/extra_info.txt"), append = T)

  line <- paste0("Elapsed time: ", time.fifference[3])
  write(line, paste0(outDir, "/extra_info.txt"), append = T)
  line <- paste0("Memory used at start: ", memory[1])
  write(line, paste0(outDir, "/extra_info.txt"), append = T)

  line <- paste0("Memory used after data cleaning: ", memory[2])
  write(line, paste0(outDir, "/extra_info.txt"), append = T)

  line <- paste0("Memory used after first consensus: ", memory[3])
  write(line, paste0(outDir, "/extra_info.txt"), append = T)

  line <- paste0("Memory used after UMI merging: ", memory[4])
  write(line, paste0(outDir, "/extra_info.txt"), append = T)

  line <- paste0("Memory used after second consensus: ", memory[5])
  write(line, paste0(outDir, "/extra_info.txt"), append = T)

  # Collapsed Fastq (no names are provided to the collapsed reads)
  file <- ShortReadQ(
    DNAStringSet(consensus_mean$read1),
    FastqQuality(consensus_mean$quality1)
    # ,
    # BStringSet(paste0(newUMIs$ID1, " ", consensus_mean$UMI))
  )
  
  message("ShortRead object generated!")

  fileSplit <- as.data.table(str_split(fastq, "\\/"))
  fileSplit <- as.data.table(str_split(fileSplit[nrow(fileSplit)], "\\."))
  output <- paste0(outDir, "/", fileSplit[1], "_corrected.fastq.gz")
  part <- fileSplit[1]
  file.create(output)
  writeFastq(file, output, mode = "a")

  message("ShortRead object saved!")

  ## For now, these lines are commenting unless we can add the names to the
  ## sequences

  # output.csv <- as.data.table(cbind(
  #   UMI = consensus_mean$UMI,
  #   UMIs = newUMIs$UMI,
  #   counts = newUMIs$Counts,
  #   read1 = consensus_mean$read1,
  #   quality1 = consensus_mean$quality1,
  #   ID1 = paste0(newUMIs$ID1, " ", consensus_mean$UMI)
  # ))
  #
  # write.table(output.csv,
  #             paste0(outDir, "/", part, "_summary_table.csv"),
  #             sep = "\t", row.names = F)

  # remove(part, file, output, fileSplit, output.csv)
}



.groupingSingle <- function(intermediate.table,
                            full,
                            quality,
                            UMIlength) {
  intermediate.table.c1 <- intermediate.table[
    which(intermediate.table$count == 1),
  ]
  intermediate.table.c2 <- intermediate.table[
    which(intermediate.table$count != 1),
  ]

  rm(intermediate.table)

  res1 <- data.table()

  if (nrow(intermediate.table.c1) > 0) {
    f1 <- full[which(full$UMI %in% intermediate.table.c1$UMI), ]
    f1 <- f1[order(f1$id), ]

    qr1 <- quality[which(quality$id %in% f1$id), ]
    qr1 <- qr1[order(qr1$id), ]
    qr1 <- as.matrix(qr1[, 1:(ncol(qr1) - 1)])
    qr1 <- qr1 + 33
    qr1 <- base::apply(qr1, 1, intToUtf8)
    qr1 <- as.character(qr1)

    result.1 <- data.table(
      UMI = f1$UMI12,
      read1 = f1$read, quality1 = qr1
    )

    colnames(result.1) <- c("UMI", "read1", "quality1")

    rm(f1, qr1)

    res1 <- result.1
  }

  if (nrow(intermediate.table.c2) > 0) {
    res2 <- plyr::ldply(
      parallel::mclapply(intermediate.table.c2$UMI, function(x) {
        r1 <- x

        # reads with specific UMI
        grouping <- full[which(full$UMI == r1), ]

        quality.1 <- quality[which(quality$id %in% grouping$id), ]

        grouping <- grouping[order(grouping$id), ]

        quality.1 <- quality.1[order(quality.1$id), ]

        # File 1
        grouping_q <- cbind(
          grouping,
          quality.1[, grep(
            pattern = "id",
            colnames(quality.1),
            invert = T
          )]
        )

        rm(grouping, quality.1)

        result1 <- .calculationsFunction(grouping_q = grouping_q)

        result.2 <- data.table(
          UMI = substr(r1, 1, UMIlength),
          read1 = result1[1, 1], quality1 = result1[1, 2]
        )

        colnames(result.2) <- c("UMI", "read1", "quality1")

        return(result.2)
      },
      mc.preschedule = FALSE,
      mc.cores = parallel::detectCores()
      ),
      data.table,
      .id = NULL
    )
  }

  result <- rbindlist(list(res1, res2))

  return(result)
}


.groupingFinalSingle <- function(newUMIs, # r1,
                                 full,
                                 quality,
                                 first_consensus,
                                 UMIlength) {
  newUMIs.1 <- newUMIs[which(str_length(newUMIs$UMI) == UMIlength), ]
  newUMIs.2 <- newUMIs[which(str_length(newUMIs$UMI) != UMIlength), ]

  rm(newUMIs)

  res1 <- data.table()

  if (nrow(newUMIs.1) > 0) {
    result.1 <- first_consensus[which(first_consensus$UMI %in% newUMIs.1$UMI), ]
    colnames(result.1) <- c("UMI", "read1", "quality1")

    res1 <- result.1
  }

  if (nrow(newUMIs.2) > 0) {
    res2 <- plyr::ldply(
      parallel::mclapply(newUMIs.2$UMI, function(x) {
        r1 <- x

        grouping <- full[str_detect(full$UMI, as.character(r1)), ]

        quality.1 <- quality[which(quality$id %in% grouping$id), ]

        grouping <- grouping[order(grouping$id), ]

        quality.1 <- quality.1[order(quality.1$id), ]

        # File 1
        grouping_q <- cbind(grouping, quality.1[, grep(
          pattern = "id", x = colnames(quality.1), invert = TRUE
        )])

        rm(grouping, quality.1)

        result1 <- .calculationsFunction(grouping_q)

        result.2 <- data.table(
          UMI = substr(r1, 1, 12),
          read1 = result1[1, 1], quality1 = result1[1, 2]
        )

        colnames(result.2) <- c("UMI", "read1", "quality1")

        res2 <- result.2
        return(res2)
      },
      mc.preschedule = FALSE,
      mc.cores = parallel::detectCores()
      ),
      data.table,
      .id = NULL
    )
  }

  result <- rbindlist(list(res1, res2))

  return(result)
}


.UMIcorrectionSingle <- function(intermediate.table,
                                 first_consensus,
                                 sequenceDistance,
                                 UMIdistance,
                                 outDir) {
  dist.calc <- 0

  uniqueUMIs <- c()
  IDs_1 <- c()
  counts <- c()

  while (nrow(intermediate.table) > 1) {
    best <- first_consensus[
      which(first_consensus$UMI == intermediate.table$UMI[1]),
    ]
    list.best <- best$UMI[1]
    list.counts <- intermediate.table$count[1]

    temp.intermediate <- intermediate.table[2:nrow(intermediate.table), ]

    base_dist <- stringdist::stringdistmatrix(
      a = best$UMI[1],
      b = temp.intermediate$UMI,
      method = "hamming"
    )[1, ]

    dist.calc <- dist.calc + length(base_dist)

    who <- which(base_dist <= UMIdistance)

    if (length(who) > 0) {
      temp.intermediate <- temp.intermediate[who, ]

      temp_read1 <- first_consensus[
        which(first_consensus$UMI %in% temp.intermediate$UMI),
      ]$read1

      dist1 <- stringdist::stringdistmatrix(
        a = best$read1,
        b = temp_read1,
        method = "hamming"
      )[1, ]

      dist.calc <- dist.calc + (length(dist1))
      who <- which(dist1 <= sequenceDistance)

      if (length(who) > 0) {
        temp.intermediate <- temp.intermediate[who, ]
        list.best <- c(list.best, temp.intermediate$UMI)
        list.counts <- c(list.counts, temp.intermediate$count)

        list.best <- paste(list.best, collapse = "|")
        list.counts <- paste(list.counts, collapse = "|")
      }
    }


    IDs_1 <- append(IDs_1, intermediate.table$ID1[1])
    counts <- append(counts, list.counts)
    uniqueUMIs <- append(uniqueUMIs, list.best)
    intermediate.table <- intermediate.table[str_detect(intermediate.table$UMI,
      as.character(list.best),
      negate = T
    ), ]
  }

  if (!is.null(intermediate.table)) {
    IDs_1 <- append(IDs_1, intermediate.table$ID1[1])
    counts <- append(counts, intermediate.table$count[1])
    uniqueUMIs <- append(uniqueUMIs, intermediate.table$UMI[1])
  }

  newUMIs <- as.data.table(cbind(
    UMI = uniqueUMIs,
    ID1 = IDs_1,
    Counts = counts
  ))

  line <- paste0("Number of Hamming distances calculated: ", dist.calc)
  write(line, paste0(outDir, "/extra_info.txt"), append = T)
  return(newUMIs)
}

.one.run.calculationsFunction <- function(one.base) {
  one.base <- one.base[, .(
    mean = mean(V2),
    count = .N,
    perc_qual = mean(V2) / 93 * 100,
    perc_count = .N / nrow(one.base) * 100
  ), by = list(V1)]

  one.base$criterion <- rowMeans(one.base[, c("perc_qual", "perc_count")])

  ### what happens if we have the same maximum criterion for different V1 ##

  one.base <- one.base[which(one.base$criterion == max(one.base$criterion)), ]

  return(data.table(
    V1 = as.character(one.base[which.max(one.base$perc_qual), ]$V1),
    mean = round(one.base[which.max(one.base$perc_qual), ]$mean)
  ))
}

.calculationsFunction <- function(grouping_q) {

  # list with all nts
  cons <- str_split(grouping_q$read, pattern = "", simplify = TRUE)
  cons <- as.list(as.data.table(cons))

  # list with all qualities
  grouping_q <- as.list(grouping_q[, 4:ncol(grouping_q)])

  # merge lists
  grouping_q <- base::Map(data.table, cons, grouping_q)

  rm(cons)

  ##
  cons_corr <- parallel::mclapply(grouping_q, function(x) 
    .one.run.calculationsFunction(one.base = x), mc.preschedule = FALSE, 
    mc.cores = detectCores())
  
  cons_corr <- rbindlist(cons_corr)

  # join again in one final sequence
  consensus <- str_c(cons_corr$V1, collapse = "")
  meanQuality <- as.numeric(cons_corr$mean) + 33
  meanQuality <- intToUtf8(meanQuality)
  result <- data.table(seq = consensus, qual = meanQuality)

  return(result)
}
