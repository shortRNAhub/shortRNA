

## longRNAs to tree
longRNA2FL <- function(file, root, cores = detectCores()) {
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
  library(data.tree)
  library(parallel)

  gff <- readGFF(file)
  suppressWarnings(gff <- data.frame(gff) %>%
    mutate_all(funs(replace_na(., ""))))
  gff$gID <- gsub("\\..*", "", gff[, "gene_id"])
  gff$tID <- gsub("\\..*", "", gff[, "transcript_id"])

  gff.sp <- split(x = gff, f = gff[, "gene_type"])

  gff.sp[[1]]$pathString <- .paste5(root, "/", gff.sp[[1]][, "gene_type"], "/", gff.sp[[1]][, "gID"], ":", gff.sp[[1]][, "gene_name"],
    "/", gff.sp[[1]][, "tID"],
    na.rm = T, sep = ""
  )

  # gff1 <- gff d <- split(gff1,rep(1:16,each=113500)) gff <- d[[1]] pathString <- mclapply(d, function(gff){ .paste5(root,
  # '/', gff[, 'gene_type'], '/', gff[, 'gID'], ':', gff[, 'gene_name'], '/', gff[, 'tID'], na.rm = T, sep = '' ) },
  # mc.preschedule = F, mc.cores = 16)

  # pathString1 <- as.character(unlist(pathString))


  node <- as.Node(data.frame(gff.sp[[1]]), na.rm = TRUE)

  gff.tree <- mclapply(gff.sp[2:length(gff.sp)], function(x) {
    x$pathString <- .paste5(x[, "gene_type"], "/", x[, "gID"], ":", x[, "gene_name"], "/", x[, "tID"], na.rm = T, sep = "")
    node <- as.Node(data.frame(x), na.rm = TRUE)
    return(node)
  }, mc.preschedule = F, mc.cores = cores)

  tmp <- lapply(gff.tree, function(x) node$AddChildNode(x))

  return(node)
}

# 
# anno <- readRDS("../../shortRNA_reports/sperm_gapp2014/shortRNA/03_tse/annotation/anno.rds")
# features <- data.frame(anno@elementMetadata)
# 
# tm <- features[features$transcript_type %in% c("tRNA"),]
# 
# x <- tm$transcript_id
# root <-  "tRNA"

## tRNA/ miRNA
tRNAmiRNAToFL <- function(features, root) {
  x <- features
  x <- sapply(strsplit(gsub(".", "-", x, fixed = TRUE), "-", fixed = TRUE), FUN = function(x) {
    if (x[1] == "mmu") {
      x <- x[-1]
    }
    n <- NULL
    if (root == "miRNA") {
      if (length(grep(pattern = "[0-9][a-z]|[0-9][0-9][a-z]", x = x[2])) > 0) {
        a <- strsplit(x[2], "")[[1]]
        a1 <- grep(pattern = "[0-9]|[0-9][0-9]", x = a, value = T)
        a2 <- grep(pattern = "[a-z]", x = a, value = T)
        n <- c(a1, paste0(a1, a2))
        x <- c(x[1], n, x[3:length(x)])
        x1 <- sapply(1:2, FUN = function(i) paste(x[1:i], collapse = "-"))
        x2 <- sapply(3:length(x), FUN = function(i) paste(x[c(1, 3:i)], collapse = "-"))
        paste(c(x1, x2), collapse = "/")
      }
    } else if (root == "tRNA") {
      x <- sapply(2:length(x), FUN = function(i) paste(x[1:i], collapse = "-"))
      paste(x, collapse = "/")
    }
  })
  ps <- paste0(root, "/", x)
  fl <- FactorList(strsplit(x, "/"))
}


## rRNAs
.rRNA2path <- function(files, root) {
  library(Biostrings)
  library(stringr)
  library(data.tree)
  files <- list.files(path = files, pattern = ".fa", full.names = T)
  files.read <- lapply(files, readDNAStringSet)

  rr <- str_extract(string = files, pattern = "[0-9]+S|[0-9].[0-9]S")
  names(files.read) <- rr

  pathString <- paste(root, names(files.read), sep = "/")
  return(as.Node(data.frame(pathString)))
}

#' Convert the annotation fasta file into tree obtained from `prepareAnnotation` function
#' @author Deepak Tanwar (tanward@ethz.ch)
#' @import Biostrings data.tree plyr
#' @seealso Biostrings data.tree plyr
#' @param fasta a link to the fasta file
#' @param root root name: longRNA or miRNA or tRNA
#' @return A `list`: a `data.frame` and a `data.tree`
#' \dontrun{
#' @examples
#' tRNA <- fastaToTree(file = "test/tRNAs.modified.fa", root = "tRNA")
#' miRNA.mirBase <- fastaToTree(file = "test/miRNAs.modified.fa", root = "miRNA")
#' rRNA <- fastaToTree(file = "test/Mus_musculus/rRNAdb/", root = "rRNA")
#' longRNA.gencode <- fastaToTree(file = "test/gencode.vM22.chr_patch_hapl_scaff.annotation.gff3.gz", root = "longRNA")
#' @export
featureToFL <- function(anno, root) {
  if (root == "longRNA") {
    return(.longRNA2path(file = file, root = root, cores = 4))
  } else if (root == "rRNA") {
    return(.rRNA2path(files = file, root = "rRNA"))
  } else {
    library(Biostrings)
    library(data.tree)
    fasta <- readDNAStringSet(file)
    features <- sapply(strsplit(names(fasta), " "), FUN = function(x) x[1])
    pathString <- tRNAmiRNAToFL(features = features, root = root)
    node <- as.Node(data.frame(features = features, pathString = pathString), na.rm = T)
    return(node)
  }
}


# Adding reads to features ----
.objToString <- function(obj) {
  deparse(substitute(obj))
}