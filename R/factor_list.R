# library(data.tree)
# # tRNAs as FL
# load("../../../shortRNA_data/db/tRNA.rda")
# ps_tRNA <- ToDataFrameTable(tRNA, "pathString")
# names(ps_tRNA) <- gsub(pattern = ".*\\/", replacement = "", x = ps_tRNA)
# 
# # Data subset
# ar_tRNA <- readRDS("ar_tRNA.rds")
# 
# 
# # Example
# ps <- ps_tRNA
# mappedFeaturesDF <- ar_tRNA
# featuresCol <- "transcript_id"
# readsCol <- "seq"

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

  if (sm > 0) warning(paste(sm, "features were removed as they were not found in the tree!"))

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

# tree_tRNA <- addReadsToTree(ps = ps_tRNA, 
#                             mappedFeaturesDF = ar_tRNA, 
#                             featuresCol = "transcript_id", 
#                             readsCol = "seq", 
#                             root = "mm10")


# library(TreeSummarizedExperiment)
# 
# m <- readRDS("reads_with_counts.rds")
# rownames(m) <- m$seq
# m <- m[,grep(pattern = "sample", x = colnames(m))]
# cd <- data.frame(Samples = colnames(m), Group = rep(c("CTRL", "MSUS"), rep = 6))
# rownames(cd) <- cd$Samples
# 
# tse <- TreeSummarizedExperiment(assay = list(counts = m), rowTree = tree_tRNA, colData = cd)
# 
# library(treeclimbR)
# res <- runDA(TSE = tse, feature_on_row = TRUE, assay = 1, 
#              option = "glm", group_column = "Group",
#              design_terms = "Group")
# out <- nodeResult(object = res, n = Inf)