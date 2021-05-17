#' mergeAtomicLists
#'
#' Pairwise merging of two atomic lists of the same length. Does not perform
#' any matching of the names, but triggers an error if names are present and
#' not in the same order.
#'
#' @param x An `AtomicList` object
#' @param y An `AtomicList` object
#'
#' @return An `AtomicList` object of the same length as `x` and `y`
#' @export
#'
#' @importFrom IRanges FactorList AtomicList
#' @importFrom S4Vectors splitAsList
#' @examples
#' f1 <- FactorList(list(LETTERS[1:3], LETTERS[1:4]))
#' f2 <- FactorList(list(LETTERS[4:6], LETTERS[7:9]))
#' mergeAtomicLists(f1, f2)
mergeAtomicLists <- function(x, y) {
  stopifnot(length(y) == length(x))
  if (is.null(names(x)) || is.null(names(y))) {
    n <- seq_along(x)
    if (!is.null(names(x))) {
      n <- names(x)
    }
    if (!is.null(names(y))) {
      n <- names(y)
    }
  } else {
    stopifnot(identical(names(x), names(y)))
    n <- names(x)
  }
  if (!is.integer(n)) {
    n <- as.factor(n)
  }
  if (is(x, "FactorList") && is(y, "FactorList")) {
    lvls <- union(levels(x)[[1]], levels(y)[[1]])
    fx <- factor(unlist(x, use.names = FALSE), lvls)
    fy <- factor(unlist(y, use.names = FALSE), lvls)
    return(splitAsList(factor(c(fx, fy), levels = seq_along(lvls), lvls), c(rep(n, lengths(x)), rep(n, lengths(y)))))
  } else {
    if (is(x, "FactorList")) {
      x <- CharacterList(x)
    }
    if (is(y, "FactorList")) {
      y <- CharacterList(y)
    }
  }
  splitAsList(c(unlist(x, use.names = FALSE), unlist(y, use.names = FALSE)), c(rep(n, lengths(x)), rep(n, lengths(y))))
}

#' longestOrderedOverlap
#'
#' Returns the longest ordered overlap between elements of an atomic list
#'
#' @param x An AtomicList
#'
#' @return An atomic vector
#' @export
#' @importFrom IRanges AtomicList
#'
#' @examples
#' il <- IntegerList(list(1:4, 1:5, c(1:3, 6)))
#' longestOrderedOverlap(il)
longestOrderedOverlap <- function(x) {
  if (!is(x, "AtomicList")) {
    if (is(x, "list")) {
      x <- AtomicList(x)
    } else {
      stop("x should be an atomic list")
    }
  }
  x1 <- x[[1]]
  i <- rep(seq_along(x), lengths(x))
  n <- length(x)
  x <- unlist(x)
  d <- which(!duplicated(i))
  j <- 0
  while (length(x) > 0 && length(d) == n && length(u <- unique(x[d])) == 1) {
    j <- j + 1
    x <- x[-d]
    i <- i[-d]
    d <- which(!duplicated(i))
  }
  x1[seq_len(j)]
}



#' fList2tree
#'
#' Converts a `FactorList` to a `phylo` object.
#'
#' @param fL A `FactorList` object
#' @param addRoot Logical; whether to add a root node
#' @param collapseSingles Logical; wether to remove intermediate nodes of single
#' childs
#'
#' @return A `phylo` object
#' @export
fList2tree <- function(fL, addRoot = TRUE, collapseSingles = FALSE, root = "ROOT") {
  fL <- relist(droplevels(unlist(fL, use.names = FALSE)), fL)
  lvls <- levels(fL[[1]])
  # temporarily get rid of levels
  x <- IntegerList(fL)
  maxLen <- max(lengths(x))
  f <- rep(seq_along(x), lengths(x))
  x <- unlist(x)
  w <- which(!duplicated(f))
  xfrom <- integer(0)
  xto <- integer(0)
  while (length(w) > 0) {
    xfrom <- c(xfrom, x[w])
    x <- x[-w]
    f <- f[-w]
    w <- which(!duplicated(f))
    xto <- c(xto, x[w])
    wdup <- which(f %in% unique(f[-w]))
    x <- x[wdup]
    f <- f[wdup]
    w <- which(!duplicated(f))
  }
  m <- matrix(c(xfrom, xto), nrow = length(xfrom))
  rm(xfrom, xto)
  m <- m[!duplicated(m), , drop = FALSE]
  if (addRoot) {
    roots <- setdiff(unique(m[, 1]), unique(m[, 2]))
    m <- rbind(cbind(rep(length(lvls) + 1L, length(roots)), roots), m)
    colnames(m) <- NULL
    lvls <- c(lvls, root)
  }
  tips <- setdiff(unique(m[, 2]), unique(m[, 1]))
  ro <- unique(c(tips, max(m), unique(m[, 1])))
  conv <- seq_along(ro)[order(ro)]
  m <- matrix(conv[as.integer(m)], ncol = 2)
  lvls <- lvls[ro]
  tips <- setdiff(unique(m[, 2]), unique(m[, 1]))
  nodes <- setdiff(unique(m[, 1]), tips)
  o <- list()
  class(o) <- "phylo"
  o$edge <- m
  o$Nnode <- length(nodes)
  o$node.label <- lvls[nodes]
  o$tip.label <- lvls[tips]
  if (collapseSingles) {
    o <- ape::collapse.singles(tree = o, root.edge = TRUE)
  }
  o
}


.addNames2fL <- function(fL) {
  nlvls <- length(levels(fL)[[1]])
  lvls <- c(levels(fL)[[1]], names(fL))
  x <- c(as.integer(unlist(fL, use.names = FALSE)), nlvls + seq_len(length(fL)))
  x <- factor(x, seq_len(nlvls + length(fL)), lvls)
  splitAsList(x, c(rep(seq_along(fL), lengths(fL)), seq_along(fL)))
}


# library(data.tree) # tRNAs as FL load('../../../shortRNA_data/db/tRNA.rda') ps_tRNA <- ToDataFrameTable(tRNA, 'pathString')
# names(ps_tRNA) <- gsub(pattern = '.*\\/', replacement = '', x = ps_tRNA) # Data subset ar_tRNA <- readRDS('ar_tRNA.rds')
# # Example ps <- ps_tRNA mappedFeaturesDF <- ar_tRNA featuresCol <- 'transcript_id' readsCol <- 'seq'

addReadsToTree <- function(ps, mappedFeaturesDF, featuresCol = "transcript_id", readsCol = "seq", ...) {
  suppressPackageStartupMessages({
    library(data.tree)
    library(plyr)
    library(future.apply)
    library(data.tree)
    library(IRanges)
    library(plyr)
    library(parallel)
  })

  # Parallel processing of apply functions
  plan(multisession)

  # pathString to factorList
  fl <- FactorList(strsplit(ps, "/"))

  # Features per sequence
  features <- strsplit(mappedFeaturesDF[, featuresCol], ";")
  featuresPerRead <- features
  names(featuresPerRead) <- mappedFeaturesDF[, readsCol]


  # Filtering out non-existing features
  featuresPerRead <- lapply(featuresPerRead, function(x) intersect(x = x, y = names(fl)))
  sm <- sum(lengths(featuresPerRead) == 0)

  if (sm > 0) {
    warning(paste(sm, "features were removed as they were not found in the tree!"))
  }

  featuresPerRead <- featuresPerRead[lengths(featuresPerRead) > 0]

  # Features per sequence factorList
  fpr_fl <- FactorList(featuresPerRead)

  # For sequences overlapping with 1 feature
  fpr_fl_s <- fpr_fl[lengths(fpr_fl) == 1]

  fpr_fl_s_o <- fl[unlist(fpr_fl_s)]
  names(fpr_fl_s_o) <- names(fpr_fl_s)
  fpr_fl_s_o <- .addNames2fL(fpr_fl_s_o)

  fpr_s_tree <- fList2tree(fL = fpr_fl_s_o, addRoot = T, collapseSingles = F)

  # For sequences overlapping with more than 1 feature
  fpr_fl_m <- fpr_fl[lengths(fpr_fl) > 1]
  fli <- IntegerList(fl)
  fpr_fl_m_o <- future_lapply(fpr_fl_m, FUN = function(x) longestOrderedOverlap(fli[x]))
  rm(fli)

  fpr_fl_m_o1 <- splitAsList(
    x = factor(as.numeric(unlist(fpr_fl_m_o)), levels = seq_along(levels(fl[[1]])), labels = levels(fl[[1]])),
    f = rep(names(fpr_fl_m_o), lengths(fpr_fl_m_o))
  )

  fpr_fl_m_o1 <- .addNames2fL(fpr_fl_m_o1)

  # Features for the tree
  fpr_tree <- c(fpr_fl_s_o, fpr_fl_m_o1)

  tree <- fList2tree(fL = fpr_tree, addRoot = T, collapseSingles = F, ...)

  plan(sequential)

  return(tree)
}

# tree_tRNA <- addReadsToTree(ps = ps_tRNA, mappedFeaturesDF = ar_tRNA, featuresCol = 'transcript_id', readsCol = 'seq', root
# = 'mm10')


# library(TreeSummarizedExperiment) m <- readRDS('reads_with_counts.rds') rownames(m) <- m$seq m <- m[,grep(pattern =
# 'sample', x = colnames(m))] cd <- data.frame(Samples = colnames(m), Group = rep(c('CTRL', 'MSUS'), rep = 6)) rownames(cd)
# <- cd$Samples tse <- TreeSummarizedExperiment(assay = list(counts = m), rowTree = tree_tRNA, colData = cd)
# library(treeclimbR) res <- runDA(TSE = tse, feature_on_row = TRUE, assay = 1, option = 'glm', group_column = 'Group',
# design_terms = 'Group') out <- nodeResult(object = res, n = Inf)
