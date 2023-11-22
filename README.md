<img src="logo/baseplot.png" width="180" height="180" />

# `shortRNA` project
This is a project for the development of short RNA analysis R package. It was initiated by Pierre-Luc in 2017.

## `shortRNA` analysis steps in `R` (a short vignette)

### Package installation
```{r, eval = FALSE}
install.packages("BiocManager")
BiocManager::install("remotes")
remotes::install_github("shortRNAhub/shortRNA")
library(shortRNA)
```

### Quality check and trimming
```{r, eval = FALSE}
qc_SE(file = fq_files, outdir = "output/", ad1 = adapter_sequence)
```

### FastQ files to sequence by counts matrix
```{r, eval = FALSE}
m <- fastq2SeqCountMatrix(files = trimmed_fastq_files)
```

### Unique sequences for alignment
```{r, eval = FALSE}
fa <- DNAStringSet(row.names(m))
names(fa) <- paste0("S", 1:length(fa))
writeXStringSet(fa, fasta_file)
```


### Obtaining databases for analysis
```{r, eval = FALSE}
db <- getDB()
```


### Annotation preparation and genome index generation for alignment
```{r, eval = FALSE}
a <- prepareAnnotation(
  ensdb = db$ensdb,
  output_dir = genomeDir,
  extra.gr = list(piRNA = db$piRNA_GR, miRNA = db$miRNA_GR),
  extra.seqs = list(rRNA = db$rRNA_fa, tRNA = db$tRNA_fa),
  resolveSplicing = NULL,
  rules = defaultAssignRules(),
  tRNAEnsembleRemove = FALSE,
  clusterMiRNA = TRUE
)
```


### Alignment
```{r, eval = FALSE}
alignShortRNA(
  fastq = "unique.fasta",
  index = "genomeDir/customGenome",
  outDir = "align", GTF = exonsBy(db$ensdb),
  GTF.featureType = "exon", GTF.attrType = "gene_id"
)
```


### Overlapping aligned data with annotations
```{r, eval = FALSE}
o <- overlapWithTx2(
  bamFile = align_file, annotation = a,
  ignoreStrand = TRUE, nbthreads = 16
)
```


### Assignment of overlapping reads to features with assignment rules
```{r, eval = FALSE}
ar <- assignReads(sources = o, rules = defaultAssignRules())
```


### Converting features to `FactorList`
```{r, eval = FALSE}
fl <- featuresAnnoToFL(a)
names(fl) <- lapply(fl, function(x) as.character(x[length(x)]))
```

### Assignment of reads to the features tree
```{r, eval=FALSE}
mappedFeaturesDF <- ar
mappedFeaturesDF$seq <- rownames(mappedFeaturesDF)
ar_tree <- addReadsToTree(
  fL = fl,
  mappedFeaturesDF = mappedFeaturesDF,
  unassigned = FALSE,
  extraTreeBranch = NULL
)
```


### Creation of TreeSummarizedExperiment object
```{r, eval = FALSE}
library(TreeSummarizedExperiment)

rt <- ar_tree
as <- list(counts = m)
cd <- DataFrame(samples = c(s1, s2), group - c("group1", "group2"))

tse <- TreeSummarizedExperiment(
  assays = as,
  rowTree = rt,
  colData = cd,
  rowData = ar[row.names(m), ],
  metadata = list(
    assignedReads = ar,
    counts = m,
    notAligned = getReadsFromBam(
      bam = align_file
    )
  )
)
```


### Differential analysis
```{r, eval = FALSE}
library(treeClimbR)
dea <- runDA(TSE = tse, filter_min_count = 20)
out <- nodeResult(object = dea, n = Inf)

cand <- getCand(
  tree = dea$tree,
  score_data = out, node_column = "node",
  p_column = "PValue", sign_column = "logFC",
  message = FALSE, threshold = 0.05
)

candB <- evalCand(
  tree = dea$tree,
  levels = cand$candidate_list,
  score_data = out, node_column = "node",
  p_column = "PValue", sign_column = "logFC",
  method = "BH", limit_rej = 0.05,
  use_pseudo_leaf = FALSE,
  message = FALSE
)

result <- topNodes(object = candB, n = Inf, p_value = 0.05)
```
