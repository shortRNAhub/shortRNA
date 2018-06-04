setClass(
  "shortRNAexp",
  slots = representation(
    seqcounts =           "matrix",
    description =         "character",
    sources =             "data.frame",
    phenoData =           "data.frame",
    composition =         "list",
    norm =                "list",
    alignStats =          "list",
    agcounts =            "matrix",
    agcounts_ambiguous =  "matrix",
    agdef =               "data.frame",
    creationDate =        "Date",
    rules =               "list",
    allsrcs =             "data.frame",
    objvers =             "numeric"
    ),
  prototype = prototype(description=NA_character_, alignStats=list(), agcounts=matrix(), agcounts_ambiguous=matrix(), 
                      agdef=data.frame(), norm=list(norm.factors=NA_real_, lib.sizes=NA_real_, cnq=NULL), allsrcs=data.frame(), 
                      composition=list(), rules=list(), creationDate=Sys.Date(), objvers=1), 
  validity=function(object){
    errors <- c()
    if(!all(colnames(object@seqcounts)==row.names(object@phenoData))) errors <- c(errors, "The colnames of `seqcounts` should correspond to row.names of `phenoData`.")
    if(!all( c("sequence","cigar","location","pos_in_feature","feature_length","src_id","src_name","src_type","overlap") %in% colnames(object@sources) )) errors <- c(errors, '`sources` should have at least the following columns: "sequence","cigar","location", "pos_in_feature", "feature_length", "src_id","src_name","src_type","overlap"')
    if(!is.numeric(object@seqcounts[1,1])) errors <- c(errors, "`seqcounts` should be a numeric matrix.")
    if(length(errors)==0) return(TRUE)
    errors
  }
)


.recastSources <- function(sources){
  sources <- as.data.frame(sources, stringsAsFactors=F)
  sources$location <- as.character(sources$location)
  sources$src_name <- as.character(sources$src_name)
  sources$overlap <- suppressWarnings(as.integer(as.numeric(as.character(sources$overlap))))
  sources$length <- suppressWarnings(as.integer(as.numeric(as.character(sources$length))))
  if("pos_in_feature" %in% colnames(sources))   sources$pos_in_feature <- suppressWarnings(as.integer(as.numeric(as.character(sources$pos_in_feature))))
  if("feature_length" %in% colnames(sources))   sources$feature_length <- suppressWarnings(as.integer(as.numeric(as.character(sources$feature_length))))
  return(sources)
}

setMethod("initialize", "shortRNAexp", function(.Object, ...) {
  o <- callNextMethod(.Object, ...)
  validObject(o)
  
  if(is.null(o@rules[["miRNA"]])) o@rules <- defaultAssignRules()
  
  sr <- o@sources[which(o@sources$sequence %in% row.names(o@seqcounts)), c("sequence","cigar","location","pos_in_feature","feature_length","src_id","src_name","src_type","overlap")]
  sr$src_name <- gsub("mmu-miR-","mmu-mir-",gsub("p(A|C|G|T)$","p",sr$src_name),fixed=T)
  sr$length <- sapply(sr$sequence, FUN=nchar)
  
  sr <- .recastSources(sr)
  
  message("
    Assigning reads to features...
")
  
  DEBUG <- FALSE
  
  if(DEBUG){
    par <- NULL
  }else{
    par <- .checkPara()
  }
  
  if(is.null(par)){
    sources <- t(sapply(unique(sr$sequence), sr=sr, rules=o@rules, DEBUG=DEBUG, FUN=function(s,sr,rules,DEBUG){ 
      a <- try(assignRead(sr[which(sr$sequence==s),,drop=F],rules),silent=T)
      if(is(a,"try-error")){
        if(DEBUG){
          print("Error trying to assign:")
          print(sr[which(sr$sequence==s),,drop=F])
        }
        stop(a[[1]])
      }
      return(a[-5])
    }))
    colnames(sources) <- c("sequence", "location", "cigar", "pos_in_feature", "src_name", "src_type", "overlap", "length", "status")
  }else{
    fns <- c("assignRead", ".disambiguateOverlappingFeatures", ".disambiguate_miRNAs", ".disambiguate_tRNA", "defaultAssignRules")
    sources <- foreach(s=unique(sr$sequence) , .export=fns, .combine=rbind) %dopar% {
      res <- as.data.frame(assignRead(sr[which(sr$sequence==s),,drop=F],o@rules)[-5])
      names(res) <- c("sequence", "location", "cigar", "pos_in_feature", "src_name", "src_type", "overlap", "length", "status")
      res
    }
    stopCluster(par)
  }
  row.names(sources) <- sources[,1]
  
  sources <- sources[,-1]
  sources <- .recastSources(sources)
  
  o@sources <- sources
  
  o@alignStats <- .getAlignStats(o, sr)
  rm(sr)
  
  message("
Aggregating sequences into features:")
  o <- aggregateSequences(o)
  
  o@composition$raw <- estimateComposition(o)
  
  return(o)
})

setGeneric("phenoData", function(object, ...){ standardGeneric("phenoData", ...) })
setMethod("phenoData","shortRNAexp", function(object){ object@phenoData })

setMethod("show", "shortRNAexp", function(object){
  message(paste0("A `shortRNAexp` object of ",nrow(object@seqcounts)," small RNA sequences (",sum(object@sources$status=="unique")," of which can be uniquely assigned) across ",ncol(object@seqcounts)," samples."))
})


setMethod("[","shortRNAexp", function(x, i){
  x@phenoData <- x@phenoData[i,,drop=F]
  x@seqcounts <- x@seqcounts[,i,drop=F]
  x@agcounts <- x@agcounts[,i,drop=F]
  x@agcounts_ambiguous <- x@agcounts_ambiguous[,i,drop=F]
  x@composition$raw$all <- x@composition$raw$all[,i,drop=F]
  x@composition$raw$abridged <- x@composition$raw$abridged[,i,drop=F]
  x@alignStats$reads <- x@alignStats$reads[,i,drop=F]
  if(all(!is.na(x@norm$norm.factors))){
    x@composition$normalized$all <- x@composition$normalized$all[,i,drop=F]
    x@composition$normalized$abridged <- x@composition$normalized$abridged[,i,drop=F]
    x@norm$norm.factors <- x@norm$norm.factors[i]
    x@norm$lib.sizes <- x@norm$lib.sizes[i]
    for(s in c("offset","glm.offset","agoffset","agglmoffset")){
      if(!is.null(x@norm[[s]])) x@norm[[s]] <- x@norm[[s]][,i,drop=F]
    }
  }
  return(x)
})

setMethod("names","shortRNAexp", function(x){
  colnames(x@seqcounts)
})

setMethod("names<-","shortRNAexp", function(x, value){
  if(length(value) != length(unique(value))) stop("Some names are in duplicate!")
  if(length(value) != ncol(x@seqcounts)) stop("The number of names given does not match the number of samples.")
  row.names(x@phenoData) <- value
  colnames(x@seqcounts) <- value
  colnames(x@agcounts) <- value
  colnames(x@agcounts_ambiguous) <- value
  colnames(x@composition$raw$all) <- value
  colnames(x@composition$raw$abridged) <- value
  colnames(o@alignStats$reads) <- value
  if(!is.null(x@composition$normalized)){
    colnames(x@composition$normalized$all) <- value
    colnames(x@composition$normalized$abridged) <- value
  }
  return(x)
})
