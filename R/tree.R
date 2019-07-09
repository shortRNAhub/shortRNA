# Convert `fasta` to `data.tree` ----

.RNAtoTree <- function(RNA) {
  library(plyr)
  n <- lapply(RNA$names, function(x) {
    a <- strsplit(x, "-")[[1]]
    if (a[1] == "mmu") a <- a[-1]
    df <- data.frame(matrix(ncol = length(a)), stringsAsFactors = F)
    for (i in 1:length(a)) {
      df[1, i] <- paste(a[1:i], collapse = "-")
    }
    if (length(grep(pattern = "\\.", x = df[, ncol(df)])) == 1) {
      b <- strsplit(df[, ncol(df)], "\\.")[[1]]
      df[, ncol(df)] <- b[1]
      df[, ncol(df) + 1] <- paste(b[1:2], collapse = ".")
    }
    colnames(df) <- c("R", paste0("C", 1:(ncol(df) - 1)))
    df$pathString <- paste(df[1, ], collapse = "/")
    return(df)
  })
  n1 <- ldply(n, data.frame)
  return(data.frame(RNA, n1, stringsAsFactors = F))
}


#' Convert the annotation fasta file into tree obtained from `prepareAnnotation` function
#' @author Deepak Tanwar (tanward@ethz.ch)
#' @import Biostrings data.tree plyr
#' @seealso Biostrings data.tree plyr
#' @param fasta a link to the fasta file
#' @return A `list`: a `data.frame` and a `data.tree`
#' \dontrun{
#' @examples
#' tRNA <- fastaToTree(fasta = "../test/tRNAs.modified.fa")
#' miRNA <- fastaToTree(fasta = "../test/miRNAs.modified.fa")
#' }
#' @export
fastaToTree <- function(fasta) {
  library(Biostrings)
  library(data.tree)
  read.file <- readDNAStringSet(fasta)
  df <- data.frame(names = names(read.file), Sequence = as.character(read.file), stringsAsFactors = F)
  df$names <- as.character(sapply(df$names, function(x) strsplit(x, "\\ ")[[1]][1]))
  df <- .RNAtoTree(RNA = df)
  tree.df <- as.Node(df)
  return(list(df = df, tree = tree.df))
}



# Combine two `data.frame` into one ----

#' Combine a list of `data.frame` into one `data.tree`. `data.frame` should have a column "pathString", or specify.
#' @author Deepak Tanwar (tanward@ethz.ch)
#' @import data.tree plyr
#' @seealso data.tree plyr
#' @param dfList a list of data frames
#' @param rootName name for the root
#' @return A `list`: a `data.frame` and a `data.tree`
#' \dontrun{
#' @examples
#' mm10 <- combineFeaturesDF(dfList = list(tRNA$df, miRNA$df), rootName = "mm10")
#' }
#' @export
combineFeaturesDF <- function(dfList, rootName = "root", ...) {
  library(plyr)
  library(data.tree)
  df <- Reduce(function(df1, df2) rbind.fill(df1, df2), dfList)
  df$pathString <- paste(rootName, df$pathString, sep = "/")
  tree.df <- as.Node(df)
  return(list(df = df, tree = tree.df))
}



# Adding reads to features ----

.objToString <- function(obj) {
  deparse(substitute(obj))
}

#' Add reads to the features in the annotation `data.tree`. Also add a separate node with "Unassigned" reads/features
#' @author Deepak Tanwar (tanward@ethz.ch)
#' @import data.tree
#' @seealso data.tree
#' @param tree a `data.tree` object
#' @param mappedFeaturesDF a `data.frame` with "Reads" and mapped "Features" columns.
#' @param featuresCol Column with features
#' @param readsCol column with reads
#' @return A `data.tree` object. Modified version of original tree input.
#' \dontrun{
#' @examples
#' # Input
#' testDF <- data.frame(
#'   Reads = c("Read1", "Read1", "Read1", "Read2", "Read3", "Read4", "Read5"),
#'   Features = c(
#'     "tRNA-Ala-XYZ-1", "tRNA-Ala-XYZ-2", "tRNA-Ala-ABC-1",
#'     "tRNA-Sel-GCG-1", "tRNA-Arg-TCG-3-1", "tRNA-Arg-TCG-3-1", "tRNA-Arg-TCG-3-1"
#'   ),
#'   stringsAsFactors = F
#' )
#'
#' # Analysis
#' testRun <- addReadsFeatures(tree = mm10$tree, mappedFeaturesDF = testDF)
#'
#' # Output
#' testTree <- FindNode(node = mm10$tree, name = "tRNA-Arg-TCG")
#' testTreeUnassigned <- FindNode(node = mm10$tree, name = "Unassigned")
#' }
#' @export
addReadsFeatures <- function(tree, mappedFeaturesDF, featuresCol = "Features", readsCol = "Reads") {
  library(data.tree)
  notFound <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = F)
  colnames(notFound) <- colnames(mappedFeaturesDF)
  for (i in 1:nrow(mappedFeaturesDF)) {
    if (is.null(FindNode(node = tree, name = mappedFeaturesDF[i, featuresCol]))) {
      message("Coud not find the feature ", mappedFeaturesDF[i, featuresCol])
      notFound[nrow(notFound) + 1, ] <- mappedFeaturesDF[i, ]
    } else {
      a <- FindNode(node = tree, name = mappedFeaturesDF[i, featuresCol])$pathString
      a <- strsplit(a, "\\/")[[1]]
      a <- a[2:length(a)]
      a <- paste(a, collapse = "`$`")
      a <- paste0("`", a, "`")
      expr <- paste0(.objToString(obj = tree), "$", a, "$AddChild(", mappedFeaturesDF[i, readsCol], ")")
      expr <- paste0(.objToString(tree), "$", a)
      eval(parse(text = expr))$AddChild(mappedFeaturesDF[i, readsCol])
    }
  }
  notFound$pathString <- paste("Unassigned", notFound[, featuresCol], notFound[, readsCol], sep = "/")
  tree$AddChildNode(child = as.Node(notFound))
  return(tree)
}