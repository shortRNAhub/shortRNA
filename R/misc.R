## Function to check for the no. of cores and for parallel processing
.checkPara <- function(ncores = NULL, maxCores = 8) {
  if (!is.null(ncores)) {
    if (ncores == 1) {
      return(NULL)
    }
  }
  if (tryCatch(require("doParallel", character.only = TRUE),
    error = function(e) FALSE
  )) {
    library(doParallel)
    if (is.null(ncores)) {
      ncores <- min(detectCores() - 1, maxCores)
    }
    if (ncores > 1) {
      cl <- makeForkCluster(nnodes = ncores, outfile = "")
      registerDoParallel(cl)
      return(cl)
    }
  }
  return(NULL)
}



#' autoLayout
#'
#' Creates a layout for a given (minimum) number of frames
#'
#' @param nb Minimum number of frames
#' @param byrow Whether to fill frames by row instead of by column
#'  (default FALSE, i.e. fill by column)
#'
#' @return The total number of frames in the layout.
#'
#' @export
autoLayout <- function(nb, byrow = FALSE) {
  nc <- ceiling(sqrt(nb))
  nr <- ceiling(nb / nc)
  layout(matrix(1:(nr * nc), nrow = nr, byrow = byrow))
  return(nr * nc)
}



## Funtion for calculating GC content in sequence
.gcContents <- function(x) {
  sapply(as.character(x), FUN = function(x) {
    x <- strsplit(x, "", fixed = T)[[1]]
    sum(x %in% c("G", "C")) / length(x)
  })
}


# finds the longest common substring
LCS <- function(x) {
  a <- sapply(as.character(x), split = "", fixed = T, FUN = strsplit)
  a <- sapply(a, mm = min(sapply(a, length)), FUN = function(x, mm) {
    x[1:mm]
  })
  w <- which(!apply(a, 1, FUN = function(x) {
    length(unique(x)) == 1
  }))[1]
  paste(a[1:w, 1], collapse = "")
}


#' capitalizeRead
#'
#' Alters a sequence read capitalization according to it's
#' alignment characteristics (cigar).
#'
#' @param read A character vector of length 1 containing the sequnce.
#' @param cigar A character vector of length 1 containing the cigar.
#' @param indels Logical; whether to also indicate insertions/deletions
#'  (default FALSE). This cannot be enablied if the read already includes
#'   any hyphen ('-').
#'
#' @return A character vector of length 1 containing the formatted sequence.
#'
#' @export

capitalizeRead <- function(read, cigar, indels = FALSE) {
  if (is.na(cigar) | !grepl("S|H|X|I", cigar)) {
    return(read)
  }
  read <- toupper(as.character(read))
  cigar <- .splitCigar(cigar) ## Please see the .splitCigar function below
  read <- strsplit(read, "", fixed = T)[[1]]
  if (any(read == "-")) {
    indels <- FALSE
    hyph <- which(read == "-")
    read <- read[-hyph]
  } else {
    hyph <- NULL
  }
  pos <- as.numeric(cigar[, 2])
  for (i in 1:nrow(cigar)) {
    x1 <- ifelse(i == 1, 1, sum(pos[1:(i - 1)]))
    if (cigar[i, 1] %in% c("H", "S", "X", "I")) {
      x <- x1:(x1 + pos[i])
      read[x] <- tolower(read[x])
    } else {
      if (indels & cigar[i, 1] == "D") {
        read <- c(read[1:(x1 - 1)], rep("-", pos[i]), read[x1:length(read)])
      }
    }
  }
  if (!is.null(hyph)) {
    .rehyphenate(read, hyph)
  } ## Please see the .rehypernate function below
  return(paste(read, collapse = ""))
}



## Function to split the cigar string.  This is used in `capitalizeRead` function
.splitCigar <- function(cigar) {
  cigar <- as.character(cigar)
  p <- gregexpr("([0-9]*[MSXIH])", cigar)[[1]]
  cigar <- strsplit(cigar, "", fixed = T)[[1]]
  p <- as.numeric(p) + attr(p, "match.length") - 1
  t(sapply(1:length(p), cigar = cigar, p = p, FUN = function(x, p, cigar) {
    c(cigar[p[x]], paste(cigar[(ifelse(x == 1, 1, p[x - 1] + 1)):as.numeric(p[x] - 1)], collapse = ""))
  }))
}


## Function to rehypernate the read This is used in `capitalizeRead` function

.rehyphenate <- function(r, h) {
  r2 <- r
  for (f in h) {
    if (length(r2) < f) {
      r2 <- c(r2[1:(f - 1)], "-")
    } else {
      if (f == 1) {
        r2 <- c("-", r2)
      } else {
        r2 <- c(r2[1:(f - 1)], "-", r2[f:(length(r2))])
      }
    }
  }
  r2
}

# use instead of `stop` to log the error and traceback
.fstop <- function(x, object = NULL) {
  futile.logger::flog.error(x)
  futile.logger::flog.trace(traceback(4))
  object <- futile.logger::ftry(object, error = function(x) NULL)
  if (!is.null(object)) {
    futile.logger::flog.trace("Offending object:", object, capture = TRUE)
  }
  stop(x, call. = FALSE)
}



.fcheck <- function(x, fatal = TRUE, object = NULL) {
  cmd <- deparse(substitute(x))
  x <- tryCatch(x,
    warning = function(w) futile.logger::flog.warn(w),
    error = function(e) {
      .fstop(paste0(
        "Error in `", cmd,
        "`: ", e$message
      ))
    }
  )
  msg <- paste0("(", cmd, ") = ", x)
  if (x || !fatal) {
    futile.logger::flog.trace(msg)
    return(x)
  }
  .fstop(msg, object = tryCatch(object, error = function(e) NULL))
}


.fstop <- function(x, object = NULL) {
  cmd <- deparse(substitute(x))
  x <- futile.logger::ftry(x, error = function(e) stop(e$message, call. = TRUE))
  if (!is.null(x) && length(x) == 1 && is.logical(x)) {
    if (x) {
      return(invisible(TRUE))
    }
    x <- paste0("(", cmd, ") == ", x)
  }
  futile.logger::flog.error(x)
  futile.logger::flog.trace(traceback(4))
  object <- futile.logger::ftry(object, error = function(e) NULL)
  if (!is.null(object)) {
    futile.logger::flog.trace("Offending object:", object, capture = TRUE)
  }
  stop(x, call. = FALSE)
}


#' Run `lapply` function in parallel, if possible.
#' @author Deepak Tanwar (tanward@ethz.ch)
#'
#' @import parallel
#'
#' @param any parameters for `lapply` or `mclapply`
#' @return A `list`.
#'
#' @export
mylapply <- function(...) {
  if (require(parallel) && .Platform$OS.type == "unix") {
    mclapply(..., mc.cores = detectCores(), mc.preschedule = F)
  } else {
    lapply(...)
  }
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

#' List files on a FTP server
#' @import RCurl stringr
#'
#' @param url A url of FTP location
#'
#' @return A list of files
#'
#' @examples
#' # Input
#' url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M18/"
#'
#' # Analysis
#' files <- listFilesFTP(url)
#'
#' # Output
#' files
listFilesFTP <- function(url) {
  library(RCurl)
  library(XML)
  # library(stringr)
  files <- getHTMLLinks(getURL(
    url = url,
    ftp.use.epsv = FALSE,
    dirlistonly = TRUE
  ))
  
  files <- paste0(url, files)
  return(files)
}
