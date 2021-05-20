#' prepareAnnotation
#'
#' @param ensdb An 'EnsDb' object, as obtained from AnnotationHub
#' @param genome A genome fasta, BSGenome or TwoBit object. If missing, will be
#' fetched using `ensdb`
#' @param output_dir The path where to save the files
#' @param extra.gr A list of GRanges or GRangesList giving the coordinates of
#' extra features. If a GRanges containing a 'type' column, exons will be used
#' and split by transcript. The objects should minimally contain a `tx_id`
#' field; if `tx_biotype` is missing and `extra.gr` is named, the names will be
#' used as biotype. `symbol` can also be provided.
#' @param extra.seqs Extra transcript sequences to include. Can either be a
#' named list of named characters (in which case the names of the list are used
#' as biotype, and the names of the character vector as tx_id), or a
#' `DNAStringSet` object with `mcols` columns 'tx_id' and 'tx_biotype' (`symbol`
#' can also be provided).
#' @param resolveSplicing The biotypes (in the `ensdb`) for which to resolve
#' splicing. If NULL, will resolve all but the low-priority (i.e. long RNA)
#' transcripts.
#' @param ... passed to `Rsubread::buildindex`
#'
#' @return Produces a Rsubread index and returns a GRanglesList of features
#' @export
#'
#' @import S4Vectors GenomicRanges Biostrings
#' @importFrom Rsubread buildIndex
prepareAnnotation <- function(ensdb, genome = NULL, output_dir = "",
                              extra.gr = list(), extra.seqs = NULL,
                              resolveSplicing = NULL, ...) {

  # Loading libraries early in case they are not available,
  # program will fail early
  library(dplyr)
  library(ensembldb)
  library(R.utils)
  library(Rsubread)

  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # If splicing is not to be resolved
  if (is.null(resolveSplicing)) {
    p <- defaultAssignRules()$highPriorityTypes
    resolveSplicing <- names(p)[p >= 0]
  }

  # If no additional sequence files are provided
  if (!is.null(extra.seqs)) {
    stopifnot(is.list(extra.seqs) || is(extra.seqs, "DNAStringSet"))

    # If additional sequence files are provided as a list
    if (is.list(extra.seqs)) {
      m <- dplyr::bind_rows(
        lapply(extra.seqs, FUN = function(x) {
          data.frame(
            row.names = names(x), tx_id = names(x),
            seq = as.character(x)
          )
        }),
        .id = "tx_biotype"
      )

      extra.seqs <- DNAStringSet(m$seq)
      m$seq <- NULL
      names(extra.seqs) <- row.names(m)
      mcols(extra.seqs) <- m
    }
    m <- mcols(extra.seqs)

    if (is.null(m$tx_id)) m$tx_id <- names(extra.seqs)
    names(extra.seqs) <- paste0("pseudoChr_", names(extra.seqs))
    stopifnot(all(c("tx_id", "tx_biotype") %in% colnames(mcols(extra.seqs))))

    if (is.null(m$symbol)) m$symbol <- m$tx_id
    gr <- GRanges(names(extra.seqs), IRanges(1L, width = nchar(extra.seqs)),
      strand = "+", tx_id = m$tx_id, tx_biotype = m$tx_biotype,
      symbol = m$symbol
    )
  }

  if (!is.null(resolveSplicing)) {
    gs <- genes(ensdb,
      filter = ~ !(tx_biotype %in% resolveSplicing),
      columns = c("gene_id", "gene_biotype", "symbol")
    )
  } else {
    gs <- genes(ensdb, columns = c("gene_id", "gene_biotype", "symbol"))
  }

  gs$tx_id <- gs$symbol
  gs$tx_biotype <- gs$gene_biotype
  gs$gene_id <- gs$gene_biotype <- NULL

  if (!is.null(resolveSplicing)) {
    tx <- exonsBy(ensdb,
      filter = ~ tx_biotype %in% resolveSplicing,
      columns = c("tx_id", "tx_biotype", "symbol")
    )
  } else {
    tx <- exonsBy(ensdb, columns = c("tx_id", "tx_biotype", "symbol"))
  }

  tx <- c(tx, split(gs, gs$tx_id))

  if (!is.null(extra.seqs)) tx <- c(tx, split(gr, gr$tx_id))

  if (length(extra.gr) > 0) {
    extra.gr <- lapply(names(extra.gr), FUN = function(x) {
      gr <- extra.gr[[x]]
      if (is(gr, "GRangesList")) {
        stopifnot(!is.null(names(gr)))
        if (!("tx_biotype" %in% colnames(gr@unlistData))) {
          gr@unlistData$tx_biotype <- x
        }
        if (!("symbol" %in% colnames(gr@unlistData))) {
          gr@unlistData$symbol <- gr@unlistData$tx_id
        }
        gr@unlistData <- gr@unlistData[, c("tx_id", "tx_biotype", "symbol")]
      } else {
        stopifnot("tx_id" %in% colnames(mcols(gr)))
        if (!("tx_biotype" %in% colnames(mcols(gr)))) mcols(gr)$tx_biotype <- x
        if (!("symbol" %in% colnames(mcols(gr)))) mcols(gr)$symbol <- x
        if (!is.null(gr$type) && any(gr$type == "exon")) gr <- gr$exon
        mcols(gr) <- mcols(gr)[, c("tx_id", "tx_biotype", "symbol")]
        gr <- split(gr, gr$tx_id)
      }
      gr
    })


    tx <- c(tx, do.call(c, extra.gr))
    anno.out <- file.path(output_dir, "features.rds")
  }

  saveRDS(tx, file = anno.out)
  message("Features saved in \n", anno.out)

  # If genome is nuot provided
  if (is.null(genome)) genome <- getSeq(ensembldb::getGenomeTwoBitFile(ensdb))
  if (is.character(genome) && length(genome) != 1) {
    stop(
      "`genome` should be the path to a fasta(.gz) file, or a ",
      "BSgenome/TwoBitFile object"
    )
  }


  if (is(genome, "TwoBitFile") || grepl(pattern = "\\.2bit", x = genome)) {
    genome <- import(genome)
  }

  isCompressed <- !is.character(genome) || grepl("\\·gz$", genome)
  genome.out <- paste0("customGenome.fasta", ifelse(isCompressed, ".gz", ""))
  genome.out <- file.path(output_dir, genome.out)


  if (is.null(extra.seqs)) { # in case of no extra seq
    if (!is.character(genome)) { # if genome is not a file
      writeXStringSet(genome, genome.out, compress = TRUE)
    } else { # if genome is a file
      file.copy(genome, genome.out)
      if (!isCompressed) R.utils::gzip(genome.out)
    }
  } else { # in case of extra seq
    if (!is.character(genome)) { # if genome is not a file
      genome <- c(genome, extra.seqs)
      writeXStringSet(genome, genome.out, compress = TRUE)
    } else { # If genome is a file
      writeXStringSet(extra.seqs,
        filepath = genome.out,
        append = TRUE
      )
      if (!isCompressed) R.utils::gzip(genome.out)
    }
  }

  rm(genome)
  message(
    "Genome including eventual extra chromosomes was saved in:\n",
    genome.out, "\nNow building the index..."
  )
  Rsubread::buildindex(basename = paste0(output_dir, "/customGenome"), 
                       reference = ifelse(test = isCompressed, 
                                          yes = genome.out,
                                          no = paste0(genome.out, ".gz")),  ...)
  return(tx)
}


ah <- AnnotationHub()
ensdb <- rev(query(ah, "GRCm38", "Ensdb"))[[1]]

tRNA <- readDNAStringSet("../../ssc/annotation/mm10-mature-tRNAs.fa")
names(tRNA) <- trimws(gsub(pattern = ".*_|\\(.*", replacement = "", x = names(tRNA)))

piRNA <- import("../extdata/piRNA_precursors_zamore_mm10.gtf")
piRNA <- piRNA[piRNA$type == "transcript"]
colnames(mcols(piRNA))[5:7] <- c("symbol", "tx_id", "tx-type")



test1 <- prepareAnnotation(
  ensdb = ensdb,
  genome = "/mnt/IM/tmp/ensdb/reference.fa",
  output_dir = "/mnt/IM/tmp/ensdb/test1",
  extra.seqs = list(tRNA = tRNA),
  extra.gr = list(piRNA = piRNA)
)


test2 <- prepareAnnotation(
  ensdb = ensdb,
  genome = "/mnt/IM/tmp/ensdb/reference.2bit",
  output_dir = "/mnt/IM/tmp/ensdb/test2",
  extra.seqs = list(tRNA = tRNA),
  extra.gr = list(piRNA = piRNA)
)

test3 <- prepareAnnotation(
  ensdb = ensdb,
  genome = "/mnt/IM/tmp/ensdb/reference.2bit",
  output_dir = "/mnt/IM/tmp/ensdb/test3",
  extra.seqs = list(tRNA = tRNA)
)

test4 <- prepareAnnotation(
  ensdb = ensdb,
  genome = "/mnt/IM/tmp/ensdb/reference.2bit",
  output_dir = "/mnt/IM/tmp/ensdb/test4",
  extra.gr = list(piRNA = piRNA)
)

test5 <- prepareAnnotation(
  ensdb = ensdb,
  genome = "/mnt/IM/tmp/ensdb/reference.2bit",
  output_dir = "/mnt/IM/tmp/ensdb/test5"
)


test6 <- prepareAnnotation(
  ensdb = ensdb,
  output_dir = "/mnt/IM/tmp/ensdb/test6"
)


test7 <- prepareAnnotation(
  ensdb = ensdb,
  output_dir = "/mnt/IM/tmp/ensdb/test7",
  extra.seqs = list(tRNA = tRNA)
)

test8 <- prepareAnnotation(
  ensdb = ensdb,
  output_dir = "/mnt/IM/tmp/ensdb/test8",
  extra.gr = list(piRNA = piRNA)
)

test9 <- prepareAnnotation(
  ensdb = ensdb,
  output_dir = "/mnt/IM/tmp/ensdb/test9",
  extra.seqs = list(tRNA = tRNA),
  extra.gr = list(piRNA = piRNA)
)