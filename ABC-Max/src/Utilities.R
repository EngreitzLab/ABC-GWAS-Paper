## Jesse Engreitz
## August 20, 2012
## R utilities

grepany <- function(patterns, x, ...) {
  sort(unlist(sapply(patterns, function(pattern) grep(pattern, x, ...))))
}

greplany <- function(patterns, v) {
  match <- rep(FALSE, length(v))
  for (pattern in patterns) {
    match <- match | grepl(pattern, v)
  }
  return(match)
}

readRefGene <- function(file) {
  x <- read.delim(file)
  colnames(x) <- c("geneName","refseq","chr","strand","start","end","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","symbol","cdsStartStat","cdsEndStat","exonFrames")
  return(x)
}

readShortBed <- function(file, chr=NULL, sort=F) {
  return(readBed(file, extra.colnames=c(), chr=chr, sort=sort))
}

readBedgraph <- function(file, chr=NULL, sort=F) {
  return(readBed(file, extra.colnames="score", chr=chr, sort=sort))
}

bed.extra.colnames <- c("name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
readBed <- function(file, extra.colnames=bed.extra.colnames, chr=NULL, sort=F, skip.chr.sorting=FALSE) {
  first.line <- readLines(file, n=1)
  if (regexpr("track", first.line) != -1) skip=1 else skip=0
  result <- read.delim(file, stringsAsFactors=T, header=F, comment.char="#", skip=skip)
  colnames(result) <- c("chr","start","end", extra.colnames)[1:ncol(result)]
  if (!is.null(chr)) result <- result[result$chr == chr,]
  if (!skip.chr.sorting) result$chr <- ordered(result$chr, levels=sortChromosomeNames(levels(result$chr)))
  if (sort) result <- result[order(result$chr, result$start, result$end),]
  return(result)    
}

readBedpe <- function(file) {
  result <- read.delim(file, stringsAsFactors=F, header=F)
  colnames(result)[1:10] <- c("chr1","start1","end1","chr2","start2","end2","name","score","strand1","strand2")
  return(result)
}

readNarrowPeak <- function(file) {
    readBed(file, extra.colnames=c("name","score","strand","signalValue","pValue","qValue","peak"))
}

sortChromosomeNames <- function(chrs) {
  x <- gsub("chr","",chrs)
  suppressWarnings(result <- c(paste("chr", sort(na.omit(as.numeric(x))), sep=''),
              paste("chr", sort(x[is.na(as.numeric(x))]), sep='')))
  return(result)
}


selectEveryNth <- function(x, n, start=1) {
  return(x[seq(start, length(x), n)])
}


GRangesFromBed <- function(bed, ...) {
  #cat("Note:  GRangesFromBed ADJUSTS start coordinate by 1 to match 0-based BED coordinates, unless start == end (different from behavior prior to 7/8/17)\n")
  bed$start <- bed$start + 1
  bed$start[bed$start > bed$end] <- bed$end[bed$start > bed$end]
  g <- GRanges(bed[,1], IRanges(bed[,2], bed[,3]), ...)
  if (ncol(bed) > 3) names(g) <- bed[,4]
  if (ncol(bed) > 4) mcols(g) <- bed[,5]
  if (ncol(bed) > 5) {
    bed[,6] <- as.character(as.matrix(bed[,6]))
    bed[bed[,6] == ".",6] <- "+"   ## GRanges doesn't like "." as strand
    strand(g) <- bed[,6]
  }
  if (ncol(bed) > 6) mcols(g) <- bed[,c(4,7:ncol(bed))]
  return(g) 
}

GRangesFromBedGraph <- function(bed) {
    #cat("Note:  GRangesFromBedgraph ADJUSTS start coordinate by 1 to match 0-based BED coordinates, unless start == end (different from behavior prior to 7/8/17)\n")
    bed$start <- bed$start + 1
    bed$start[bed$start > bed$end] <- bed$end[bed$start > bed$end]  
    g <- GRanges(bed[,1], IRanges(bed[,2], bed[,3]))
    mcols(g) <- bed[,4]
    return(g)
}

GRangesToBed <- function(gr) {
  cat("Note:  GRangesToBed ADJUSTS start coordinate by 1 to match 0-based BED coordinates (same as behavior prior to 7/8/17)\n") 
  to.keep <- 3
  df <- data.frame(chr=seqnames(gr),
                   start=start(gr)-1,
                   end=end(gr),
                   name="",
                   score=0,
                   strand='+')
  if (!is.null(gr$name)) {  
    df$name <- gr$name  
    to.keep <- 4  
  } 
  if (!is.null(gr$score)) {
    df$score <- gr$score
    to.keep <- 5
  }
  if (!is.null(strand(gr))) {
    df$strand <- as.vector(strand(gr))
    to.keep <- 6
  }
  return(df[,1:to.keep])
}


addTSSToBED <- function(bed) {
  bed$tss <- bed$start
  bed$tss[bed$strand == "-"] <- bed$end[bed$strand == "-"]
  return(bed)
}

refactor <- function(factor1, factor2) {
  factor(as.character(as.matrix(factor1)), levels=levels(factor2))
}

unfactor <- function(df) {
  for (i in 1:ncol(df)) {
    if (is.factor(df[,i])) {
      df[,i] <- as.character(as.matrix(df[,i]))
    }
  }
  return(df)
}

write.tab <- function(d, file, row.names=F, col.names=T, sep='\t', quote=F, ...) {
  write.table(d, file, row.names=row.names, col.names=col.names, sep=sep, quote=quote, ...)
}

writeBed <- function(d, file, ...) {
  tmp <- options("scipen")[[1]]
  options(scipen=999)
  write.tab(d, file, col.names=F, ...)
  options(scipen=tmp)
}

writeBed3 <- function(d, file, ...) {
  writeBed(d[,c("chr","start","end")], file, ...)
}

fixPythonLogicalFields <- function(df) {
  for (i in 1:ncol(df)) {
    if (all(df[,i] %in% c("True","False",""))) {
      df[,i] <- as.logical(as.character(as.matrix(df[,i])))
    }
  }
  return(df)
}
