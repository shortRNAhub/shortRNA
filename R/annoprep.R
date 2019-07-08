#' getAnnotationFiles
#'
#' Download (and untar) the files necessary to build a short RNA annotation, and returns a named list indicating the different feature files. Is called by `prepareAnnotation`.
#'
#' @param org An organism assembly to use. Currently supported are: mm10, hg38.
#' @param destination The destination folder for the files (defaults to current working directory)
#'
#' @return A list.
#'
#' @export
getAnnotationFiles <- function(org = "mm10", destination = getwd()) {
  if (org == "mm10") {
    ff <- list(
      c("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.tRNAs.gtf.gz", 
        file.path(destination, "gencode.tRNAs.gtf.gz")),
      c("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/gencode.vM18.chr_patch_hapl_scaff.annotation.gtf.gz", 
        file.path(destination, "gencode.features.gtf.gz")),
      c("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.tRNAs.gtf.gz",
        file.path(destination, "gencode.tRNAs.gtf.gz")),
      c("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz", 
        file.path(destination, "mirbase.fa.gz")),
      c("ftp://mirbase.org/pub/mirbase/CURRENT/genomes/mmu.gff3", 
        file.path(destination, "mirbase.gtf")),
      c("http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz", 
        file.path(destination, "gtRNAdb.tar.gz"))
    )
    srcs <- sapply(ff, FUN = function(x) {
      download.file(x[[1]], ifelse(length(x) > 1, x[[2]], NULL))
      return(x[[1]])
    })
    untar(file.path(destination, "gtRNAdb.tar.gz"))
    library(Biostrings)
    fa <- readRNAStringSet(file.path(destination, "mirbase.fa.gz"), "fasta")
    fa <- fa[grep("^mmu-", names(fa))]
    writeXStringSet(fa, file.path(destination, "mm10.mature.miRNAs.fa"), format = "fasta")
    unlink(file.path(destination, "mature.miRNA.fa.gz"))
    return(list(
      features.gtf = file.path(destination, "gencode.features.gtf.gz"),
      mirbase.fa = file.path(destination, "mm10.mature.miRNAs.fa"),
      mirbase.gtf = file.path(destination, "mirbase.gtf"),
      # pi_precursors.gtf = file.path(path.package("shortRNA"), "extdata", "mm10.piRNA.precursors.gtf"),
      pi_precursors.gtf = file.path("../inst/extdata/mm10.piRNA.precursors.gtf"),
      gtRNAdb.fa = file.path(destination, "mm10-tRNAs.fa"),
      gtRNAdb.bed = file.path(destination, "mm10-tRNAs.bed"),
      tRNA.gtf = file.path(destination, "gencode.tRNAs.gtf.gz"),
      equivalences = NA_character_,
      srcs = srcs
    ))
  }
  if (org == "hg38") {
    ff <- list(
      c("http://www.smallrnagroup.uni-mainz.de/piRNAclusterDB/homo_sapiens/proTRAC_normal_ovary_generic/BED_files.zip", file.path(destination, "piRNAs.clusters.zip")),
      c("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.tRNAs.gtf.gz", file.path(destination, "gencode.tRNAs.gtf.gz")),
      c("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.chr_patch_hapl_scaff.annotation.gtf.gz", file.path(destination, "gencode.features.gtf.gz")),
      c("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz", file.path(destination, "mirbase.fa.gz")),
      c("ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3", file.path(destination, "mirbase.gtf")),
      c("http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz", file.path(destination, "gtRNAdb.tar.gz"))
    )

    srcs <- sapply(ff, FUN = function(x) {
      download.file(x[[1]], ifelse(length(x) > 1, x[[2]], NULL))
      return(x[[1]])
    })

    unzip(file.path(destination, "piRNAs.clusters.zip"), exdir = destination)
    e <- read.delim(file.path(destination, "homo_sapiens", "proTRAC_normal_ovary_generic", "piRNA_clusters.BED"), header = F, skip = 1, stringsAsFactors = F)
    e$name <- sapply(as.character(e$V4), FUN = function(x) {
      x <- strsplit(x, "|", fixed = T)[[1]]
      paste0("locus", x[length(x)])
    })
    e$score <- 0
    e$strand <- "*"
    e$strand[grep("mono:minus", e$V4)] <- "-"
    e$strand[grep("mono:plus", e$V4)] <- "+"
    write.table(e[, -4], file.path(destination, "piRNA.clusters.bed"), col.names = F, row.names = F, sep = "\t", quote = F)
    unlink(file.path(destination, "homo_sapiens/proTRAC_normal_ovary_generic"), T)

    untar(file.path(destination, "gtRNAdb.tar.gz"))
    library(Biostrings)
    fa <- readRNAStringSet(file.path(destination, "mirbase.fa.gz"), "fasta")
    fa <- fa[grep("^hsa-", names(fa))]
    writeXStringSet(fa, file.path(destination, "hg38.mature.miRNAs.fa"), format = "fasta")
    unlink(file.path(destination, "mature.miRNA.fa.gz"))
    return(list(
      features.gtf = file.path(destination, "gencode.features.gtf.gz"),
      mirbase.fa = file.path(destination, "hg38.mature.miRNAs.fa"),
      mirbase.gtf = file.path(destination, "mirbase.gtf"),
      pi_precursors = file.path(destination, "piRNA.clusters.bed"),
      gtRNAdb.fa = file.path(destination, "hg38-tRNAs.fa"),
      gtRNAdb.bed = file.path(destination, "hg38-tRNAs.bed"),
      tRNA.gtf = file.path(destination, "gencode.tRNAs.gtf.gz"),
      equivalences = NA_character_,
      srcs = srcs
    ))
  }
  stop("Unknown organism!")
}

#' prepareAnnotation
#'
#' Prepares an annotation GRanges and modified miRNA/tRNA sequences for a given organism. Either `filelist` or `species` should be provided.
#'
#' @param species The species for automatic fetching and creation of the annotation. See `?getAnnotationFiles` for the supported species/assemblies.
#' @param filelist A list of named feature files, as produced by 'getAnnotationFiles', indicating the necessary files. The following slots are required: 'features.gtf','features.gtf','mirbase.gtf','gtRNAdb.fa','gtRNAdb.bed','tRNA.gtf'; the 'pi_precursors' is optional.
#' @param destination The destination folder for the files (defaults to current working directory)
#' @param description An optional character vector containing a description of the annotation (default empty)
#'
#' @return A GRanges object, and produces modified fasta files.
#'
#' @export
prepareAnnotation <- function(species = NULL, filelist = NULL, destination = getwd(), description = "") {
  if ((is.null(species) && is.null(filelist)) ||
    (!is.null(species) && !is.null(filelist))
  ) {
    stop("Exactly one of `species` or `filelist` must be given.")
  }
  if (!is.null(species)) {
    filelist <- getAnnotationFiles(species, destination = destination)
  }
  ff <- c("features.gtf", "features.gtf", "mirbase.gtf", "gtRNAdb.fa", "gtRNAdb.bed", "tRNA.gtf")
  if (!is.list(filelist) || !all(ff %in% names(filelist))) stop(paste("`filelist` should be a list with the following elements:", paste(ff, collapse = ", ")))
  library(Biostrings)
  library(rtracklayer)
  library(GenomicRanges)
  library(tools)

  # miRNAs
  mirs <- .buildIsomirs(filelist$mirbase.fa)
  writeXStringSet(mirs, file.path(destination, "miRNAs.modified.fa"), format = "fasta")
  mn <- sapply(names(mirs), FUN = function(x) {
    strsplit(x, " ", fixed = T)[[1]][[1]]
  })
  mn2 <- sapply(mn, FUN = function(x) {
    strsplit(x, ".", fixed = T)[[1]][[1]]
  })
  mn3 <- sapply(mn2, FUN = function(x) {
    x <- strsplit(x, "-", fixed = T)[[1]]
    paste(x[-length(x)], collapse = "-")
  })
  le <- length(mn)
  mirs <- GRanges(mn, IRanges(rep(1, le), width(mirs)), rep("+", le), transcript_id = mn2, gene_id = mn3, transcript_type = "miRNA")

  mirbase.gtf <- .prepare.mirbase.gtf(filelist$mirbase.gtf)
  mirs <- suppressWarnings(c(mirbase.gtf, mirs))

  # tRNAs
  trnas <- .preparetRNAsequences(filelist$gtRNAdb.fa, filelist$gtRNAdb.bed)
  writeXStringSet(trnas, file.path(destination, "tRNAs.modified.fa"), format = "fasta")
  le <- length(trnas)
  mn <- sapply(names(trnas), FUN = function(x) {
    x <- strsplit(x, " ", fixed = T)[[1]][[1]]
    x <- strsplit(x, "_", fixed = T)[[1]]
    x[length(x)]
  })
  trnas <- GRanges(names(trnas), IRanges(rep(1, le), width(trnas)), rep("+", le), transcript_id = names(trnas), gene_id = mn, transcript_type = "tRNA")
  tr <- import.bed(filelist$gtRNAdb.bed)
  tr$transcript_id <- tr$name
  tr$gene_id <- sapply(tr$name, FUN = function(x) {
    x <- strsplit(x, "-", fixed = T)[[1]]
    paste(x[1:3], collapse = "-")
  })
  tr$transcript_type <- "tRNA"
  g <- import.gff(filelist$tRNA.gtf)
  ov <- findOverlaps(g, tr)
  g$transcript_type <- "putative_tRNA"
  g$transcript_type[which(g$gene_type == "Pseudo_tRNA")] <- "pseudo_tRNA"
  g$gene_id <- g$transcript_name
  trnas <- suppressWarnings(c(tr[, c("transcript_id", "gene_id", "transcript_type")], trnas, g[-unique(ov@from), c("transcript_id", "gene_id", "transcript_type")]))

  # piRNAs
  if ("pi_precursors" %in% names(filelist) && !is.na(filelist$pi_precursors)) {
    pirnas <- import(filelist$pi_precursors)
    if ("type" %in% names(pirnas@elementMetadata)) {
      pirnas <- pirnas[which(pirnas$type == "transcript"), c("transcript_id", "gene_id", "transcript_type")]
    } else {
      pirnas$transcript_id <- pirnas$name
      pirnas$gene_id <- pirnas$name
      pirnas$transcript_type <- "piRNA_precursor"
      pirnas <- pirnas[, c("transcript_id", "gene_id", "transcript_type")]
    }
    trnas <- suppressWarnings(c(trnas, pirnas))
  }

  # genes
  g <- import.gff(filelist$features.gtf)
  g <- g[which(g$type %in% c("transcript", "exon")), ]
  g <- g[!duplicated(g[, c("gene_id", "transcript_id")]), ]
  twe <- unique(as.character(g$transcript_id[which(g$type == "exon")]))
  g$transcript_type[which(g$type == "transcript" && g$transcript_id %in% twe)] <- "precursor"
  g <- g[, c("transcript_id", "gene_id", "transcript_type")]

  g1 <- g[which(g$transcript_type != "miRNA"), ]
  g2 <- g[which(g$transcript_type == "miRNA"), ]
  ov <- findOverlaps(g2, mirbase.gtf)
  g2 <- g2[-unique(ov@from), ]

  if ("equivalences" %in% names(filelist) && !is.na(filelist$equivalences)) {
    gid <- sapply(as.character(g$gene_id), FUN = function(x) {
      x <- strsplit(x, ".", fixed = T)[[1]]
      if (length(x) == 1) {
        return(x)
      }
      paste(x[-length(x)], collapse = ".")
    })
    eqv <- read.delim(filelist$equivalences, header = F, row.names = 1)
    w <- which(gid %in% row.names(eqv))
    g$gene_id[w] <- eqv[gid[w], 1]
  }
  g <- suppressWarnings(c(mirs, trnas, g1, g2))

  if ("srcs" %in% names(filelist)) {
    srcs <- filelist$srcs
  } else {
    srcs <- NA
  }
  attr(g, "description") <- description
  attr(g, "srcs") <- srcs
  attr(g, "md5") <- md5sum(as.character(filelist[which(names(filelist) != "srcs" & !is.na(filelist))]))
  message(paste0("Saved modified sequences in:
", file.path(destination, "tRNAs.modified.fa"), "
", file.path(destination, "miRNAs.modified.fa")))
  g
}

.prepare.mirbase.gtf <- function(g) {
  if (is.character(g)) g <- import.gff(g)
  g$transcript_id <- g$Name
  idm <- g$Name
  names(idm) <- g$ID
  g$gene_id <- idm[g$Derives_from]
  w <- which(is.na(g$gene_id))
  g$gene_id[w] <- g$Name[w]
  g$transcript_type <- g$type
  g[, c("transcript_id", "gene_id", "transcript_type")]
}

# where mature.fa is either a RNAStringSet or a path to a fasta file of RNA sequences
.buildIsomirs <- function(mature.fa) {
  library(Biostrings)
  if (is.character(mature.fa)) mature.fa <- readRNAStringSet(mature.fa)
  fa <- DNAStringSet(complement(mature.fa))
  a <- DNAStringSet(unlist(lapply(fa, FUN = function(x) {
    sapply(c("A", "C", "G", "T"), y = x, FUN = function(y, x) {
      xscat(y, x)
    })
  })))
  names(a) <- paste(sapply(names(a), FUN = function(x) {
    strsplit(x, " ", fixed = T)[[1]][[1]]
  }), rep(c("A", "C", "G", "T"), length(a) / 4), sep = ".")
  c(fa, a)
}

.preparetRNAsequences <- function(cDNA, bed = NULL) {
  library(Biostrings)
  if (is.character(cDNA)) cDNA <- readDNAStringSet(cDNA)
  if (is.null(bed)) {
    message("No tRNA bed file provided - assuming that the fasta is stranded and spliced.")
  } else {
    # a bed file is provided, so we will splice the sequences and invert if on negative strand
    if (is.character(bed)) bed <- import.bed(bed)
    names(cDNA) <- sapply(names(cDNA), FUN = function(x) {
      x <- strsplit(x, " ", fixed = T)[[1]][[1]]
      x <- strsplit(x, "_", fixed = T)[[1]]
      x[length(x)]
    })
    names(bed) <- bed$name
    bed <- bed[names(cDNA), ]
    # splicing
    for (i in which(sapply(bed$blocks, length) > 1)) {
      cDNA[[i]] <- DNAString(paste(apply(as.data.frame(bed$blocks[[i]]), 1, fa = as.character(cDNA[[i]]), FUN = function(x, fa) {
        substr(fa, x[[1]], x[[2]])
      }), collapse = ""))
    }
    # inverting sequence
    for (i in which(as.character(strand(bed)) == "-")) {
      cDNA[[i]] <- reverse(cDNA[[i]])
    }
  }
  # we add the CCA tail
  x <- xscat(cDNA, "CCA")
  names(x) <- names(cDNA)
  x
}