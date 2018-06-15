#' getAnnotationFiles
#'
#' Download (and untar) the files necessary to build a short RNA annotation, and returns a named list indicating the different feature files.
#'
#' @param org An organism assembly to use. Currently supported are: mm10.
#'
#' @return A list.
#'
#' @export
getAnnotationFiles <- function(org="mm10"){
    if(org=="mm10"){
        ff <- list( c("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.tRNAs.gtf.gz","gencode.tRNAs.gtf.gz"),
                    c("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M14/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf.gz","gencode.features.gtf.gz"),
                    c("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M17/gencode.vM17.tRNAs.gtf.gz","gencode.tRNAs.gtf.gz"),
                    c("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz", "mirbase.fa.gz"),
                    c("ftp://mirbase.org/pub/mirbase/CURRENT/genomes/mmu.gff3","mirbase.gtf"),
                    c("http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz","gtRNAdb.tar.gz")
                )
        srcs <- sapply(ff,FUN=function(x){ download.file(x[[1]],ifelse(length(x)>1,x[[2]],NULL)); return(x[[1]]) })
        untar("gtRNAdb.tar.gz")
        library(Biostrings)
        fa <- readDNAStringSet("mirbase.fa.gz","fasta")
        fa <- fa[grep("^mmu-",names(fa))]
        writeXStringSet(fa,"mm10.mature.miRNAs.fa", format="fasta")
        unlink("mature.miRNA.fa.gz")
        return( list(features.gtf="gencode.features.gtf.gz",
                    mirbase.fa="mm10.mature.miRNAs.fa",
                    mirbase.gtf="mirbase.gtf",
                    pi_precursors.gtf=paste0(path.package("shortRNA"),"/extdata/mm10.piRNA.precursors.gtf"),
                    gtRNAdb.fa="mm10-tRNAs.fa",
                    gtRNAdb.bed="mm10-tRNAs.bed",
                    tRNA.gtf="gencode.tRNAs.gtf.gz",
                    equivalences=NA_character_,
                    srcs=srcs)
                )
    }
    stop("Unknown organism!")
}

#' prepareAnnotation
#'
#' Prepares an annotation GRanges and modified miRNA/tRNA sequences for a given organism. See `getAnnotationFiles` for obtaining the necessary files.
#'
#' @param filelist A list of named feature files, as produced by 'getAnnotationFiles', indicating the necessary files. The following slots are required: 'features.gtf','features.gtf','mirbase.gtf','gtRNAdb.fa','gtRNAdb.bed','tRNA.gtf'; the 'pi_precursors.gtf' is optional.
#'
#' @return A GRanges object, and produces modified fasta files.
#'
#' @export
prepareAnnotation <- function(filelist){
    ff <- c("features.gtf","features.gtf","mirbase.gtf","gtRNAdb.fa","gtRNAdb.bed","tRNA.gtf")
    if(!is.list(filelist) || !all(ff %in% names(filelist))) stop(paste("`filelist` should be a list with the following elements:",paste(ff,collapse=", ")))
    library(Biostrings)
    library(rtracklayer)
    library(GenomicRanges)
    library(tools)
    
    # miRNAs
    mirs <- .buildIsomirs(filelist$mirbase.fa)
    writeXStringSet(mirs,"miRNAs.modified.fa", format="fasta")
    mn <- sapply(names(mirs),FUN=function(x){ strsplit(x," ",fixed=T)[[1]][[1]] })
    mn2 <- sapply(mn, FUN=function(x){ strsplit(x,".",fixed=T)[[1]][[1]] })
    mn3 <- sapply(mn2, FUN=function(x){ x <- strsplit(x,"-",fixed=T)[[1]]; paste(x[-length(x)],collapse="-") })
    le <- length(mn)
    mirs <- GRanges(mn, IRanges(rep(1,le), width(mirs)), rep("+",le), transcript_id=mn2, gene_id=mn3, transcript_type="miRNA")
    
    mirbase.gtf <- .prepare.mirbase.gtf(filelist$mirbase.gtf)
    mirs <- suppressWarnings(c(mirbase.gtf,mirs))
    
    # tRNAs
    trnas <- .preparetRNAsequences(filelist$gtRNAdb.fa, filelist$gtRNAdb.bed)
    writeXStringSet(trnas,"tRNAs.modified.fa", format="fasta")
    le <- length(trnas)
    mn <- sapply(names(trnas), FUN=function(x){ x <- strsplit(x," ",fixed=T)[[1]][[1]]; x <- strsplit(x,"_",fixed=T)[[1]]; x[length(x)] })
    trnas <- GRanges(names(trnas), IRanges(rep(1,le), width(trnas)), rep("+",le), transcript_id=names(trnas), gene_id=mn, transcript_type="tRNA")
    tr <- import.bed(filelist$gtRNAdb.bed)
    tr$transcript_id <- tr$name
    tr$gene_id <- sapply(tr$name, FUN=function(x){ x <- strsplit(x,"-",fixed=T)[[1]]; paste(x[1:3],collapse="-") })
    tr$transcript_type <- "tRNA"
    g <- import.gff(filelist$tRNA.gtf)
    ov <- findOverlaps(g,tr)
    g$transcript_type <- "putative_tRNA"
    g$transcript_type[which(g$gene_type=="Pseudo_tRNA")] <- "pseudo_tRNA"         
    g$gene_id <- g$transcript_name
    trnas <- suppressWarnings(c(tr[,c("transcript_id","gene_id","transcript_type")],trnas,g[-unique(ov@from),c("transcript_id","gene_id","transcript_type")]))
    
    # piRNAs
    if("pi_precursors.gtf" %in% names(filelist) && !is.na(filelist$pi_precursors.gtf)){
        pirnas <- import.gff(filelist$pi_precursors.gtf)
        pirnas <- pirnas[which(pirnas$type=="transcript"),c("transcript_id","gene_id","transcript_type")]
        trnas <- suppressWarnings(c(trnas,pirnas))
    }
    
    # genes
    g <- import.gff(filelist$features.gtf)
    g <- g[which(g$type %in% c("transcript","exon")),]
    g <- g[!duplicated(g[,c("gene_id","transcript_id")]),]
    twe <- unique(as.character(g$transcript_id[which(g$type=="exon")]))
    g$transcript_type[which(g$type=="transcript" && g$transcript_id %in% twe)] <- "precursor"
    g <- g[,c("transcript_id","gene_id","transcript_type")]
    
    g1 <- g[which(g$transcript_type != "miRNA"),]
    g2 <- g[which(g$transcript_type == "miRNA"),]
    ov <- findOverlaps(g2,mirbase.gtf)
    g2 <- g2[-unique(ov@from),]
    
    if("equivalences" %in% names(filelist) && !is.na(filelist$equivalences)){    
        gid <- sapply(as.character(g$gene_id),FUN=function(x){ x <- strsplit(x,".",fixed=T)[[1]]; if(length(x)==1) return(x); paste(x[-length(x)],collapse=".") })
        eqv <- read.delim(filelist$equivalences,header=F,row.names=1)
        w <- which(gid %in% row.names(eqv))
        g$gene_id[w] <- eqv[gid[w],1]
    }
    g <- suppressWarnings(c(mirs, trnas, g1, g2))
    
    if("srcs" %in% names(filelist)){
        srcs <- filelist$srcs
    }else{
        srcs <- NA
    }
    attr(g,"info") <- list( srcs=srcs,
                            md5=md5sum(as.character(filelist[which(names(filelist)!="srcs" & !is.na(filelist))])))
                            
    message("Saved modified sequences in tRNAs.modified.fa and miRNAs.modified.fa")
    g
}

.prepare.mirbase.gtf <- function(g){
    if(is.character(g))     g <- import.gff(g)
    g$transcript_id <- g$Name
    idm <- g$Name
    names(idm) <- g$ID
    g$gene_id <- idm[g$Derives_from]
    w <- which(is.na(g$gene_id))
    g$gene_id[w] <- g$Name[w]
    g$transcript_type <- g$type
    g[,c("transcript_id","gene_id","transcript_type")]
}

# where mature.fa is either a RNAStringSet or a path to a fasta file of RNA sequences
.buildIsomirs <- function(mature.fa){
    library(Biostrings)
    if(is.character(mature.fa)) mature.fa <- readRNAStringSet(mature.fa)
    fa <- DNAStringSet(complement(mature.fa))
    a <- DNAStringSet(unlist(lapply(fa,FUN=function(x){ sapply(c("A","C","G","T"),y=x,FUN=function(y,x){ xscat(y,x)})})))
    names(a) <- paste(sapply(names(a), FUN=function(x){ strsplit(x," ",fixed=T)[[1]][[1]] }), rep(c("A","C","G","T"),length(a)/4), sep=".")
    c(fa,a)
}

.preparetRNAsequences <- function(cDNA, bed=NULL){
    library(Biostrings)
    if(is.character(cDNA)) cDNA <- readDNAStringSet(cDNA)
    if(is.null(bed)){
        message("No tRNA bed file provided - assuming that the fasta is stranded and spliced.")
    }else{
        # a bed file is provided, so we will splice the sequences and invert if on negative strand
        if(is.character(bed)) bed <- import.bed(bed)
        names(cDNA) <- sapply(names(cDNA), FUN=function(x){ x <- strsplit(x," ",fixed=T)[[1]][[1]]; x <- strsplit(x,"_",fixed=T)[[1]]; x[length(x)] })
        names(bed) <- bed$name
        bed <- bed[names(cDNA),]
        # splicing
        for(i in which(sapply(bed$blocks,length)>1)){
            cDNA[[i]] <- paste(apply(as.data.frame(bed$blocks[[i]]),1,fa=as.character(cDNA[[i]]),FUN=function(x,fa){ substr(fa,x[[1]],x[[2]])}),collapse="")
        }
        # inverting sequence
        for(i in which(strand(bed)=="-")){
            cDNA[[i]] <- reverse(cDNA[[i]])
        }
    }
    # we add the CCA tail
    x <- xscat(cDNA,"CCA")
    names(x) <- names(cDNA)
    x
}
