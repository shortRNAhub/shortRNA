#' collapseFastq
#'
#' Extract unique reads and their counts from (adapter-trimmed) fastq files, and save the contents in .seqcounts files.
#'
#' @param files A character vector of path to fastq (or fastq.gz) files.
#' @param minLength The minimum length for a read to be considered (default 15).
#' @param maxLength The maximum length for a read to be considered (default 49).
#' @param shell Shell executable, default 'bash'.

#'
#' @export
collapseFastq <- function(files, minLength=15, maxLength=49, shell="bash"){
    p <- paste0(path.package("smallRNA"),"/collapse.sh")
    f2 <- sapply(files, FUN=function(x){ paste0(gsub("\\.fastq$","",gsub("\\.gz$","",x)),".seqcounts") })
    cmds <- paste(shell,p,files,f2,minLength,maxLength)
    for(f in cmds){
        print(f)
        system(f)
    }
}

#' collapsed2countMatrix
#'
#' Creates a count matrix from a list of .seqcounts files
#'
#' @param seqcounts.files A character vector of path to .seqcounts files (i.e. files produced by `collapseFastq`).
#' @param output.file Path to the file where the output matrix will be saved.
#' @param shell Shell executable, default 'bash'.

#'
#' @export
collapsed2countMatrix <- function(seqcounts.files, output.file, shell="bash"){
    cmd <- paste0(shell," ",path.package("smallRNA"),"/collapsed2countMatrix.sh ",paste(seqcounts.files,collapse=" ")," > ",output.file)
    print(cmd)
    system(cmd)
}

#' smallRNAexp_align
#'
#' A wrapper for the 2-steps alignment method.
#'
#' @param fasta Path to a fasta file containing the unique reads.
#' @param outputfile Path to the output bam file.
#' @param bowtie1index Path to the base of the bowtie1 index.
#' @param starindex Path to the STAR index folder.
#' @param bowtie1 Executable for bowtie1 (default 'bowtie').
#' @param star Executable for STAR (default 'STAR').
#' @param samtools Executable for samtools (default 'samtools').
#' @param m Maximum alignments for a read to be considered (default 1000).
#' @param nthreads Number of threads for alignment (default 4).
#'
#' @export
smallRNAexp_align <- function(fasta, outputfile, bowtie1index, starindex, bowtie1="bowtie", star="STAR", samtools="samtools", m=1000, nthreads=4){
    cmd <- paste0(bowtie1,' -p ',nthreads,' -v 0 -S -a --best --strata -m ',m,' -f --un ',outputfile,'.unmapped.fasta ',bowtie1index,' ',fasta,' | ',samtools,' view -bh > ', outputfile, '.unsorted.bam')
    print(cmd)
    system(cmd)
    td <- tempdir()
    cmd <- paste0(samtools,' sort -T ',td,' -@ ', nthreads,' -m 2G ', outputfile,'.unsorted.bam > ', outputfile,'.align1.bam; rm ',outputfile,'.unsorted.bam')
    print(cmd)
    system(cmd)
    cmd <- paste0(star, ' --genomeDir ', starindex, ' --runThreadN ', nthreads, ' --readFilesIn ',outputfile,'.unmapped.fasta --alignIntronMax 1 --outFilterMultimapNmax ',m,' --outSAMattributes NH HI NM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ',td,'/align2 --outSAMprimaryFlag AllBestScore')
    print(cmd)
    system(cmd)
    cmd <- paste0(samtools, ' merge -f ', outputfile,' ', outputfile,'.align1.bam ', td,'align2*bam
rm ',outputfile,'.align1.bam; rm -R ',td,'
',samtools, ' index ', outputfile)
    print(cmd)
    system(cmd)
}

#' smallRNAexp_parseBam
#'
#' Extracts reads' possible sources from an alignment file.
#'
#' @param bam Path to the sorted bam file.
#' @param elements Path to the genomic elements bed file.
#' @param outputfile Path for the output file.
#' @param shell Shell executable, default 'bash'.
#' @param samtools Executable for samtools (default 'samtools').
#' @param bedtools Executable for bedtools (default 'bedtools').
#'
#' @export
smallRNAexp_parseBam <- function(bam, elements, outputfile=NULL, shell="bash", samtools="samtools", bedtools="bedtools"){
    if(is.null(outputfile)){
        f <- tempfile("srcs")
    }else{
        f <- outputfile
    }
    p <- paste0(path.package("smallRNA"),"/parseBam.sh")
    cmd <- paste0(shell,p,bam,elements,f,samtools,bedtools)
    system(cmd)
    if(is.null(outputfile)){
        ret <- read.delim(f, header=F, stringsAsFactors=F)
        unlink(f)
        return(ret)
    }
}
