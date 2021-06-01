# Convert `fasta` to `data.tree` ----

# https://stackoverflow.com/questions/43803949/create-and-print-a-product-hierarchy-tree-without-na-from-data-frame-in-r-with
.paste5 <- function(..., sep = " ", collapse = NULL, na.rm = FALSE) {
  if (na.rm == FALSE) {
    paste(..., sep = sep, collapse = collapse)
  } else
    if (na.rm == TRUE) {
      paste.na <- function(x, sep) {
        x <- gsub("^\\s+|\\s+$", "", x)
        ret <- paste(na.omit(x), collapse = sep)
        is.na(ret) <- ret == ""
        return(ret)
      }
      df <- data.frame(..., stringsAsFactors = FALSE)
      ret <- apply(df, 1, FUN = function(x) paste.na(x, sep))
      
      if (is.null(collapse)) {
        ret
      } else {
        paste.na(ret, sep = collapse)
      }
    }
}

## longRNAs to tree
.longRNA2path <- function(file, root, cores = detectCores() - 2) {
  library(rtracklayer)
  library(dplyr)
  library(tidyr)
  library(data.tree)
  library(parallel)
  
  gff <- readGFF(file)
  suppressWarnings(gff <- data.frame(gff) %>% mutate_all(funs(replace_na(., ""))))
  gff$gID <- gsub("\\..*", "", gff[, "gene_id"])
  gff$tID <- gsub("\\..*", "", gff[, "transcript_id"])
  
  gff.sp <- split(x = gff, f = gff[, "gene_type"])
  
  gff.sp[[1]]$pathString <- .paste5(root, "/", gff.sp[[1]][, "gene_type"], "/", gff.sp[[1]][, "gID"], ":",
                                    gff.sp[[1]][, "gene_name"], "/", gff.sp[[1]][, "tID"],
                                    na.rm = T, sep = ""
  )
  
  node <- as.Node(data.frame(gff.sp[[1]]), na.rm = TRUE)
  
  gff.tree <- mclapply(gff.sp[2:length(gff.sp)], function(x) {
    x$pathString <- .paste5(x[, "gene_type"], "/", x[, "gID"], ":",
                            x[, "gene_name"], "/", x[, "tID"],
                            na.rm = T, sep = ""
    )
    node <- as.Node(data.frame(x), na.rm = TRUE)
    return(node)
  }, mc.preschedule = F, mc.cores = cores)
  
  tmp <- lapply(gff.tree, function(x) node$AddChildNode(x))
  
  return(node)
}


## tRNA/ miRNA
.name2path <- function(features, root) {
  x <- features
  x <- sapply(strsplit(gsub(".", "-", x, fixed = TRUE), "-", fixed = TRUE), FUN = function(x) {
    if (x[1] == "mmu") x <- x[-1]
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
  paste0(root, "/", x)
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
#' }
#' @export
fastaToTree <- function(file, root) {
  if (root == "longRNA") {
    return(.longRNA2path(file = file, root = root, cores = 4))
  } else if(root == "rRNA"){
    return(.rRNA2path(files = file, root = "rRNA"))
  } else{
    library(Biostrings)
    library(data.tree)
    fasta <- readDNAStringSet(file)
    features <- sapply(strsplit(names(fasta), " "), FUN = function(x) x[1])
    .name2path(features = features, root = root) -> pathString
    node <- as.Node(data.frame(features = features, pathString = pathString), na.rm = T)
    return(node)
  }
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
#' Reads = c("Read1", "Read1", "Read1", "Read2", "Read3", "Read4", "Read5"),
#' Features = c(
#'  "tRNA-Ala-AGC-10", "tRNA-Ala-AGC-2", "tRNA-Ala-AGC-3",
#'  "tRNA-Ala-AGC-1C", "tRNA-Arg-AGC-1D", "tRNA-Ala-AGC-1E", "tRNA-Ala-AGC-3D"
#' ),
#' stringsAsFactors = F
#')
#'
#' data("tRNA")
#'
#' # Analysis
#' testRun <- addReadsFeatures(tree = tRNA, mappedFeaturesDF = testDF)
#'
#' # Output
#' testTree <- FindNode(node = tRNA, name = "tRNA-Ala-AGC")
#' testTreeUnassigned <- FindNode(node = mm10$tree, name = "Unassigned")
#' }
#' @export
addReadsFeatures <- function(tree,
                             mappedFeaturesDF,
                             featuresCol = "Features",
                             readsCol = "Reads") {
  # Libraries
  library(data.tree)
  library(plyr)
  library(future.apply)
  
  # Parallel processing of apply functions
  plan(multisession)
  
  features <- strsplit(mappedFeaturesDF[, featuresCol], ";")
  # featuresPerRead <- split(features, mappedFeaturesDF[, readsCol])
  featuresPerRead <- features
  names(featuresPerRead) <- mappedFeaturesDF[, readsCol]
  
  featuresPerRead.single <- featuresPerRead[lengths(featuresPerRead) == 1]
  # featuresPerRead.single$test <- "trf"
  
  featuresPerRead.notFound <- compact(lapply(featuresPerRead.single, function(x) .not_found(feature = x, tree)))
  
  featuresPerRead.single <- featuresPerRead.single[!names(featuresPerRead.single) %in% names(featuresPerRead.notFound)]
  
  featuresPerRead.multi <- featuresPerRead[lengths(featuresPerRead) > 1]
  # featuresPerRead.multi$TESTSEQ <- c("tRNA1", "tRNA2")
  
  find_parent <- .find_parent(featuresList = featuresPerRead.multi, tree = tree)
  
  # Found ones
  reads <- future_lapply(find_parent, function(x) x[[1]])
  reads <- compact(reads)
  reads <- c(reads, featuresPerRead.single)
  
  if (length(reads) > 0) {
    # Find node and add sequence
    tmp <- imap(reads, function(feature, seq) FindNode(node = tree, name = feature)$AddChild(seq))
  }
  
  # Not found ones
  notFound_single <- ldply(featuresPerRead.notFound)
  colnames(notFound_single) <- c("seq", "feature")
  
  notFound <- future_lapply(find_parent, function(x) x[[2]])
  notFound <- ldply(compact(notFound), .id = NULL)
  notFound <- rbind(notFound, notFound_single)
  
  if (nrow(notFound) > 0) {
    notFound$pathString <- paste("Unassigned", notFound[, "feature"], notFound[, "seq"], sep = "/")
    tree$AddChildNode(child = as.Node(notFound))
  }
  
  # Return the new tree with added sequences
  return(tree)
}



# Finding nodes for all the sequences
.find_parent <- function(featuresList, tree) {
  library(purrr)
  imap(featuresList, function(feature, name) {
    f <- feature
    r <- name
    
    # A dataframe for storing not found reads
    notFound <- data.frame(matrix(ncol = 2, nrow = 0), stringsAsFactors = FALSE)
    colnames(notFound) <- c("seq", "feature")
    
    p <- list()
    for (i in 1:length(f)) { # For each feature
      if (!is.null(FindNode(node = tree, name = r))) {
        message(paste("Read exist!", r, "is already assigned!"))
        next
      } else { # If read is not assigned to the feature
        if (is.null(FindNode(node = tree, name = f[i]))) { # If feature do not exists
          message("Coud not find the feature ", f[i])
          notFound[nrow(notFound) + 1, ] <- c(r, f[i])
        } else { # If feature exists
          a <- FindNode(node = tree, name = f[i])$pathString
          a <- strsplit(a, "\\/")[[1]]
          p[[i]] <- a
        }
      }
    }
    
    parent <- vector()
    
    if (length(p) == 0) { # If read do not exist
      parent <- NULL
    } else if (length(p) == 1) { # in case of 1 parent
      ca <- p[[1]]
      ca <- ca[max(na.omit(match(ca, p[[1]])))]
      parent <- ca
    } else { # In case of multiple parents
      ca <- p[[1]]
      for (i in 2:length(p)) {
        if (i != length(p)) {
          ca <- ca[na.omit(match(ca, p[[i]]))]
        } else {
          ca <- ca[max(na.omit(match(ca, p[[i]])))]
        }
      }
      parent <- ca
    }
    
    # Return features that were found and not found
    return(list(parent = parent, notAssigned = notFound))
  })
}



.not_found <- function(feature, tree){
  if (is.null(FindNode(node = tree, name = feature))) {
    return(feature)
  }
}



# Save tree children as separate data objects ----
.treeToData <- function(tree, cores = detectCores() - 2) {
  library(parallel)
  library(data.tree)
  library(usethis)
  n <- names(tree$children)
  tmp <- mclapply(n, function(x) {
    assign(x = x, value = tree[[x]])
    p <- paste0("use_data(", x, ",overwrite = TRUE)")
    eval(parse(text = p))
  }, mc.preschedule = FALSE, mc.cores = cores)
}

# .treeToData(tree = longRNA.gencode, cores = 4)

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
combineFeaturesDF <- function(dfList, rootName = "root", ...) {
  library(plyr)
  library(data.tree)
  df <- Reduce(function(df1, df2) rbind.fill(df1, df2), dfList)
  df$pathString <- paste(rootName, df$pathString, sep = "/")
  tree.df <- as.Node(df)
  return(list(df = df, tree = tree.df))
}