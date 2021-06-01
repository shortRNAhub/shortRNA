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
    if (!is.null(names(x))) n <- names(x)
    if (!is.null(names(y))) n <- names(y)
  } else {
    stopifnot(identical(names(x), names(y)))
    n <- names(x)
  }
  if (!is.integer(n)) n <- as.factor(n)
  if (is(x, "FactorList") && is(y, "FactorList")) {
    lvls <- union(levels(x)[[1]], levels(y)[[1]])
    fx <- factor(unlist(x, use.names = FALSE), lvls)
    fy <- factor(unlist(y, use.names = FALSE), lvls)
    return(splitAsList(
      c(fx, fy),
      c(rep(n, lengths(x)), rep(n, lengths(y)))
    ))
  }
  splitAsList(
    c(unlist(x, use.names = FALSE), unlist(y, use.names = FALSE)),
    c(rep(n, lengths(x)), rep(n, lengths(y)))
  )
}



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
  while (length(unique(vapply(tmp,
    FUN.VALUE = character(1),
    FUN = function(x) x[i + 1]
  ))) == 1) {
    i <- i + 1
  }
  paste(tmp[[1]][seq_len(i)], collapse = delim)
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