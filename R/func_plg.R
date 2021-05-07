#' longestCommonString
#'
#' Finds the longest common part of a set of character strings
#'
#' @param x A character vector
#' @param delim The delimiter (by default, splits all letters)
#'
#' @return A character vector of length 1
#'
#' @examples
#' a <- c("B1/B2/B3", "B1/B2/B3/B4", "B1/B2/B5")
#' longestCommonString(a, "/")
longestCommonString <- function(x, delim = "") {
  if (length(x) == 1) {
    return(x)
  }
  tmp <- strsplit(as.character(x), delim, fixed = TRUE)
  if (any(lengths(tmp) == 1)) {
    return("")
  }
  i <- 0
  while (length(unique(vapply(tmp, FUN.VALUE = character(1), FUN = function(x) x[i + 1]))) == 1) {
    i <- i + 1
  }
  paste(tmp[[1]][seq_len(i)], collapse = delim)
}


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
