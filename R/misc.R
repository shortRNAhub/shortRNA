.checkPara <- function(ncores=NULL, maxCores=8){
	if(!is.null(ncores)){
		if(ncores==1)	return(NULL)
	}
	if(tryCatch(require("doParallel", character.only=TRUE),error = function(e) FALSE)){
		library(doParallel)
		if(is.null(ncores))	ncores <- min(detectCores()-1,maxCores)
		if(ncores > 1){
			cl <- makeForkCluster(nnodes=ncores, outfile = "")
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
#' @param byrow Whether to fill frames by row instead of by column (default FALSE, i.e. fill by column)
#'
#' @return The total number of frames in the layout.
#'
#' @export
autoLayout <- function(nb, byrow = FALSE){
    nc <- ceiling(sqrt(nb))
    nr <- ceiling(nb/nc)
    layout(matrix(1:(nr * nc), nrow = nr, byrow = byrow))
    return(nr*nc)
}


.gcContents <- function(x){
    sapply(as.character(x),FUN=function(x){
        x <- strsplit(x,"",fixed=T)[[1]]
        sum(x %in% c("G","C"))/length(x)
    })
}

.splitCigar <- function(cigar){
    cigar <- as.character(cigar)
    p <- gregexpr("([0-9]*[MSXIH])",cigar)[[1]]
    cigar <- strsplit(cigar,"",fixed=T)[[1]]
    p <- as.numeric(p)+attr(p,"match.length")-1
    t(sapply(1:length(p),cigar=cigar,p=p,FUN=function(x,p,cigar){ c(cigar[p[x]], paste(cigar[(ifelse(x==1,1,p[x-1]+1)):as.numeric(p[x]-1)],collapse="")) }))
}


#' capitalizeRead
#'
#' Alters a sequence read capitalization according to it's alignment characteristics (cigar).
#'
#' @param read A character vector of length 1 containing the sequnce.
#' @param cigar A character vector of length 1 containing the cigar.
#' @param indels Logical; whether to also indicate insertions/deletions (default FALSE). This cannot be enablied if the read already includes any hyphen ('-').
#'
#' @return A character vector of length 1 containing the formatted sequence.
#'
#' @export
capitalizeRead <- function(read, cigar, indels=FALSE){
    if(is.na(cigar) | !grepl("S|H|X|I",cigar)) return(read)
    read <- toupper(as.character(read))
    cigar <- .splitCigar(cigar)
    read <- strsplit(read,"",fixed=T)[[1]]
    if(any(read=="-")){
        indels <- FALSE
        hyph <- which(read=="-")
        read <- read[-hyph]
    }else{
        hyph <- NULL
    }
    pos <- as.numeric(cigar[,2])
    for(i in 1:nrow(cigar)){
        x1 <- ifelse(i==1,1,sum(pos[1:(i-1)]))
        if(cigar[i,1] %in% c("H","S","X","I")){
            x <- x1:(x1+pos[i])
            read[x] <- tolower(read[x])
        }else{
            if(indels & cigar[i,1]=="D"){
                read <- c(read[1:(x1-1)],rep("-",pos[i]),read[x1:length(read)])
            }
        }
    }
    if(!is.null(hyph)) .rehyphenate(read, hyph)
    return(paste(read,collapse=""))
}

.rehyphenate <- function(r, h){
    r2 <- r
    for(f in h){
        if(length(r2)<f){
            r2 <- c(r2[1:(f-1)],"-")
        }else{
            if(f==1){
                r2 <- c("-",r2)
            }else{
                r2 <- c(r2[1:(f-1)],"-",r2[f:(length(r2))])
            }
        }
    }
    r2
}
