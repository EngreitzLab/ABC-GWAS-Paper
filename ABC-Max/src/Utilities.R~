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


getBedBedIntersection <- function(bed1, bed2, buffer=0) {
  bed1$chr <- as.character(bed1$chr)
  bed2$chr <- as.character(bed2$chr)
  result <- c()
  for (i in 1:nrow(bed2)) {
    result <- rbind(result, getBedIntersection(bed2[i,1], bed2[i,2], bed2[i,3], bed1, buffer))
  }
  return(result)
}

getBedIntersection <- function(chr, start, end, bed, buffer=0, getMask=FALSE, getIndices=FALSE) {
  mask <- (bed[,1] == chr) & (bed[,2] < end + buffer) & (bed[,3] > start - buffer)
  if (getMask) return(mask)
  intersectors <- which(mask)
  if (getIndices) return(intersectors)
  return(bed[intersectors,])
}


parseRegionName <- function(region) {
  region <- gsub(',','',region)
  if (length(region) == 1) {
    split1 <- strsplit(region, ":")[[1]]
    split2 <- strsplit(split1[2], "-")[[1]]
    result <- list(chr=split1[1], start=as.numeric(split2[1]), end=as.numeric(split2[2]))
    return(result)
  } else {
    result <- data.frame(t(sapply(as.character(region), parseRegionName)))
    result$chr <- unlist(result$chr)
    result$start <- as.numeric(unlist(result$start))
    result$end <- as.numeric(unlist(result$end))
    return(result)
  }
}

parseUCSC <- function(ucsc) {
  return(parseRegionName(ucsc))
}

bedToUCSC <- function(bed) {
  return(paste0(bed$chr,":",bed$start,"-",bed$end))
}


plotCumulativeDistribution <- function(v1, v2, name1, name2, ...) {
  p.value <- wilcox.test(v1, v2)$p.value
  plot(ecdf(v1), col='red', do.points=FALSE, ...)
  plot(ecdf(v2), add=TRUE, do.points=FALSE)
  legend("right", c(name1, name2), col=c("red","black"), pch=19)
  legend("bottomright", paste("Wilcoxon P-Value =", format(p.value, digits=4, scientific=T)), bty='n')
}

plotCDFs <- function(df, data.col, by.col=NULL, colors=NULL, ylab="Cumulative Frequency", legend.loc="right", legend.labels=NULL, ...) {
  ## df = data frame with 2 columns:  data.col and by.col
  ## OR:  df = data frame with multiple data.cols

  if (!is.null(by.col)) {
    if (length(data.col) > 1) stop("Cannot use multiple data.col when specifying by.col")
    labels <- sort(unique(as.character(as.matrix(df[,by.col]))))
    plot(ecdf(subset(df, get(by.col) == labels[1])[,data.col]), col=colors[1], do.points=FALSE, ylab=ylab, ...)
    for (i in 2:length(labels)) {
      plot(ecdf(subset(df, get(by.col) == labels[i])[,data.col]), col=colors[i], do.points=FALSE, add=TRUE)
    }
  } else {
    labels <- data.col
    plot(ecdf(df[,labels[1]]), col=colors[1], do.points=FALSE, ylab=ylab, ...)
    for (i in 2:length(labels)) {
      plot(ecdf(df[,labels[i]]), col=colors[i], do.points=FALSE, add=TRUE)
    }
  }
  if (is.null(legend.labels)) legend.labels <- labels
  legend(legend.loc, legend.labels, col=colors, lty=1)
}

getGenesWithinDistance <- function(lncRNA, refGene=NA, distance=1000000, fullyContained=TRUE, column="symbol") {
  if (is.na(refGene)) refGene <- readRefGene("/seq/lincRNA/data/mm9/UCSC/mm9/Annotation/Genes/refGene.txt")
  lncRNA.chr <- as.character(as.matrix(lncRNA$chr))

  if (!fullyContained) genes <- subset(refGene, as.character(as.matrix(chr)) == lncRNA.chr & (start <= lncRNA$end + distance & end >= lncRNA$start - distance))
  else genes <- subset(refGene, (as.character(as.matrix(chr)) == lncRNA.chr) & (start > lncRNA$start - distance & end < lncRNA$end + distance))
  return(as.character(as.matrix(genes[,column])))
}

splitEllipsesArguments <- function(prefixes, ...) {
  all.args <- list(...)
  result <- list()
  all.i <- c()
  for (prefix in prefixes) {
    i <- grep(paste("^",prefix,sep=''), names(all.args))
    result[[prefix]] <- all.args[i]
    all.i <- c(all.i, i)
  }
  result[["remaining"]] <- all.args[-1*all.i]
  return(result)
}

plotScaledOverlay <- function(x, y, cols=NULL, types=NULL, col='black', type='p', ...) {
  ## Plot an arbitrary number of traces (e.g. for a chromosome correlation) on top of one another
  ## ncol(y), labels, cols should all be the same length

  y.scaled <- apply(y, 2, scale)
  scaled.range <- range(y.scaled, na.rm=T)
  
  for (i in ncol(y):1) {
    if (!is.null(cols)) {
      col <- cols[i]
    }
    if (!is.null(types)) {
      type <- types[i]
    }
    
    if (i == ncol(y)) {
      plot(x, y.scaled[,i], ylim=scaled.range, col=col, type=type, yaxt='n', xaxt='n', ylab='', xlab='', ...)
    } else {
      points(x, y.scaled[,i], col=col, type=type, ...)
    }
  }
}


plotChromosomeCorrelationPanel <- function(x, y1, y2, y1.lab, y2.lab, col.y1='blue', col.y2='red', type.y1='p', type.y2='p', legend=T, x.axis=T, y.axis=T, cex=0.3, ...) {
  y1.scaled <- scale(y1)
  y2.scaled <- scale(y2)

  scaled.range <- range(c(y1.scaled,y2.scaled), na.rm=T)

  if (type.y1 == 'h' & type.y2 == 'p') {
    tmp <- y2
    y2 <- y1
    y1 <- tmp
  }
  plot(x, y2.scaled, ylim=scaled.range, col=col.y2, yaxt='n', ylab='', pch=19, cex=cex, xaxt='n', las=1, type=type.y2, ...)
  points(x, y1.scaled, col=col.y1, pch=19, cex=cex, type=type.y1)
  if (type.y1 == 'h' & type.y2 == 'p') {
    tmp <- y2
    y2 <- y1
    y1 <- tmp
  }
  
  if (legend) {
    legend("topleft", c(paste("Pearson's R =", format(cor(y1.scaled, y2.scaled, method="pearson", use="pairwise.complete.obs"), digits=2)),
                        paste("Spearman's R = ", format(cor(y1.scaled, y2.scaled, method="spearman", use="pairwise.complete.obs"), digits=2))),
           bty='n', cex=0.6)
    legend("topright", c(y1.lab, y2.lab), col=c(col.y1,col.y2), pch=19, bty='n', cex=0.6)
  }

  if (x.axis) axis(1, tck=0.02)
  if (y.axis) {
    y1.ticks <- pretty(scaled.range * sd(y1, na.rm=T) + mean(y1, na.rm=T))
    y2.ticks <- pretty(scaled.range * sd(y2, na.rm=T) + mean(y2, na.rm=T))
    
    axis(2, col.axis=col.y1, at=(y1.ticks - mean(y1, na.rm=T))/sd(y1, na.rm=T), labels=y1.ticks, tck=0.02, las=1)
    axis(4, col.axis=col.y2, at=(y2.ticks - mean(y2, na.rm=T))/sd(y2, na.rm=T), labels=y2.ticks, tck=0.02, las=1)
  }
}


plotChromosomeCorrelation <- function(x, y1, y2, y1.lab, y2.lab, smooth=1, main=paste(y1.lab,'vs',y2.lab), cex=0.7, scatter=TRUE, type='p') {
  if (scatter) layout(matrix(c(1,2), 1, 2), widths=c(1,2))
  par(cex=cex, cex.lab=cex, cex.axis=cex, cex.main=cex, mgp=c(1.5,0.3,0), oma=rep(0,4), mar=c(3,3,3,2))

  if (scatter) {
    plot(y1, y2, xlab=y1.lab, ylab=y2.lab)
    abline(lsfit(y1,y2), col='red', lty=2)
  }
  
  y1.smoothed <- filter(y1, rep(1/smooth, smooth))
  y2.smoothed <- filter(y2, rep(1/smooth, smooth))

  x.mb <- x / 1000000

  plotChromosomeCorrelationPanel(x.mb, y1.smoothed, y2.smoothed, y1.lab, y2.lab, type.y1=type, type.y2=type, xlab="Chromosome Coordinate (Mb)")
  
  #y1.scaled <- scale(y1.smoothed)
  #y2.scaled <- scale(y2.smoothed)

  
  #scaled.range <- range(c(y1.scaled,y2.scaled),na.rm=T)
  #plot(x.mb, y2.scaled, main=main, ylim=scaled.range,
  #     col='red', xlab="Chromosome Coordinate (Mb)", yaxt='n', ylab='', pch=19, cex=0.3, xaxt='n', las=1, type=type)

  #if (type == 'p') points(x.mb, y1.scaled, col='blue', pch=19, cex=0.3)
  #else if (type == 'l') lines(x.mb, y1.scaled, col='blue', pch=19, cex=0.3)

  #legend("topleft", c(paste("Pearson's R =", format(cor(y1.scaled, y2.scaled, method="pearson", use="pairwise.complete.obs"), digits=2)),
  #                    paste("Spearman's R = ", format(cor(y1.scaled, y2.scaled, method="spearman", use="pairwise.complete.obs"), digits=2))),
  #                    bty='n', cex=0.6)
  #legend("topright", c(y1.lab, y2.lab), col=c("blue","red"), pch=19, bty='n', cex=0.6)

  #axis(1, tck=0.02)

  #y1.ticks <- pretty(scaled.range * sd(y1.smoothed, na.rm=T) + mean(y1.smoothed, na.rm=T))
  #y2.ticks <- pretty(scaled.range * sd(y2.smoothed, na.rm=T) + mean(y2.smoothed, na.rm=T))
  
  #axis(2, col.axis='blue', at=(y1.ticks - mean(y1.smoothed, na.rm=T))/sd(y1.smoothed, na.rm=T), labels=y1.ticks, tck=0.02, las=1)
  #axis(4, col.axis='red', at=(y2.ticks - mean(y2.smoothed, na.rm=T))/sd(y2.smoothed, na.rm=T), labels=y2.ticks, tck=0.02, las=1)
}

odds <- function(x) {
  x[seq(1,length(x),2)]
}

plotChromosomeCorrelations <- function(regions, x, y, file, smooth=1) {
  windows <- data.frame(t(sapply(regions, parseRegionName)))
  midpoints <- unlist(windows$start) + (unlist(windows$end) - unlist(windows$start))/2

  pdf(file=file, onefile=T, width=8, height=3)

  for (j in 1:ncol(x)) {
    y1.lab <- colnames(x)[j]

    for (i in 1:ncol(y)) {
      y2.lab <- colnames(y)[i]

      tryCatch({
        plotChromosomeCorrelation(midpoints, x[,j], y[,i], y1.lab, y2.lab, smooth)
      }, error = function(e) {
        cat(paste("Skipping", y2.lab, ":", e, "\n"))
      })
    }
  }

  dev.off()
  
}


## Use this for loading a correlation matrix from CollectAnnotationCorrelations
loadCorrelationMatrix <- function(file, header=T) {
  x <- read.delim(file, header=header, stringsAsFactors=F)
  ucsc <- data.frame(t(sapply(x$X, parseUCSC)))
  x$chr <- unlist(ucsc$chr)
  x$start <- unlist(ucsc$start)
  x$end <- unlist(ucsc$end)
  x$mid <- with(x, (start + end) / 2)
  return(x)
}

## Mask out regions that are close to a given point
getSortedCorrMatrixMask <- function(matrix, point, flank=3000000) {
  return((matrix$end <= point - flank) | (matrix$start >= point + flank))
}


## Use this for loading the correlation results table from CollectAnnotationCorrelations.jar
loadCorrelationResultsTable <- function(file, remove.base="GenomeAnnotations.") {
  result <- read.delim(file, stringsAsFactors=F, header=F)
  colnames(result) <- c("analysis","feature","name","type.y1","type.y2","statistic","method","value","p","total.features","region.features")
  result$feature <- gsub(remove.base, "", result$feature)
  result <- result[order(result$value, decreasing=T),]
  return(result)
}



loadScriptureOutput <- function(file) {
  result <- read.delim(file, header=F, stringsAsFactors=F)
  colnames(result)[1:6] <- c("chr","start","end","refseq","RPKM","strand")
  colnames(result)[13] <- "p"
  result$expressed <- (result$p < 0.001)
  return(result)
}

asRPKM <- function(count, total, windowSize) {
  return(count/ (total / 1000000.0) / (windowSize / 1000.0))
}

loadSlideAndCount <- function(file) {
  curr <- read.delim(file, header=F, stringsAsFactor=F)
  colnames(curr) <- c("chr","start","end","blank","score","refseq","ratio","targetCount","controlCount","targetTotalReads","controlTotalReads","targetNormalized","controlNormalized")
  curr <- subset(curr, chr == "chrX")
  curr$mid <- (curr$start + curr$end)/2
  return(curr)
}


loadChromosomeSizes <- function(genome="mm9") {
  result <- read.delim("/seq/lincRNA/data/mm9/sizes", header=F)
  sizes <- result[,2]
  names(sizes) <- result[,1]
  return(sizes)
}

loadChromosomeSizesAsBed <- function(genome="mm9") {
  sizes <- loadChromosomeSizes()
  return(data.frame(chr=names(sizes), start=1, end=sizes))
}


bedTools.2in<-function(functionstring="bedIntersect",bed1,bed2,opt.string="") {
  ##create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  
  ##write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  
  ## create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  
  res=read.table(out,header=F)
  unlink(a.file);unlink(b.file);unlink(out)
  return(res)
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

expandBed <- function(bed.df, upstream, downstream) {
  if (!("strand" %in% colnames(bed.df)) | !all(bed.df$strand %in% c("-","+"))) {
    cat("Processing as unstranded\n")
    bed.df$start <- bed.df$start - upstream
    bed.df$end <- bed.df$end + downstream
  } else {
    cat("Processing as stranded\n")
    ## STRANDED
    pos <- bed.df$strand == "+"
    neg <- !pos

    bed.df$start[pos] <- bed.df$start[pos] - upstream
    bed.df$end[pos] <- bed.df$end[pos] + downstream

    bed.df$start[neg] <- bed.df$start[neg] - downstream
    bed.df$end[neg] <- bed.df$end[neg] + upstream
  }
  return(bed.df)
}

sortChromosomeNames <- function(chrs) {
  x <- gsub("chr","",chrs)
  suppressWarnings(result <- c(paste("chr", sort(na.omit(as.numeric(x))), sep=''),
              paste("chr", sort(x[is.na(as.numeric(x))]), sep='')))
  return(result)
}

readSlideAndCount <- function(file, chr=NULL, sort=F) {
  return(readBed(file, chr=chr, sort=sort, extra.colnames=c(bed.extra.colnames, "ratio", "ratio.log2", "count.target", "rpkm.target", "regionTotal.target", "total.target",
                         "count.control", "rpkm.control", "regionTotal.control", "total.control", "p.enriched", "p.depleted")))
}

readCountReads <- function(file, chr=NULL, sort=F) {
  return(readBed(file, chr=chr, sort=sort, extra.colnames=c(bed.extra.colnames, "count", "rpkm", "regionTotal", "total", "length")))
}

getBedMidpoints <- function(bed) {
  return(with(bed, start + (start + end)/2))
}

readExonIntronCounts <- function(exonFile, intronFile, min.coverage=0) {
  exons <- readBed(exonFile, extra.colnames=c(bed.extra.colnames, "count","fpkm","region","total"))
  introns <- readBed(intronFile, extra.colnames=c(bed.extra.colnames, "count","fpkm","region","total"))
  exons$gene <- sapply(strsplit(as.matrix(exons$name), "_exon"), "[", 1)
  introns$gene <- sapply(strsplit(as.matrix(introns$name), "_intron"), "[", 1)
  genes <- intersect(exons$gene, introns$gene)
  fpkm <- data.frame(exons=calculateFpkmFromBed(exons, min.coverage)[genes], introns=calculateFpkmFromBed(introns, min.coverage)[genes])
  fpkm$ratio <- log2((fpkm$exons + 0.1) / (fpkm$introns + 0.1))
  return(fpkm)
}

calculateFpkmFromBed <- function(count.bed, min.coverage=0) {
  with(count.bed, tapply(1:nrow(count.bed), gene, function(i) {
    coverage <- sum(count[i])
    if (coverage >= min.coverage)
      coverage / (sum(end[i] - start[i]) / 1000) / (total[1] / 1000000)
    else
      NA
  }))
}

loadInterpolatedBedgraph <- function(bedgraph.file, chr=NULL) {
  ## Return a bedgraph file with an extra column representing whether that
  ## window was an interpolated value or not
  
  interpolated.file <- sub(".bedgraph", ".interpolated.bedgraph", bedgraph.file)

  bedgraph <- readBedgraph(bedgraph.file, chr)
  bedgraph$chr <- as.factor(bedgraph$chr)
  interpolated <- readBedgraph(interpolated.file, chr)
  interpolated$chr <- as.factor(interpolated$chr)
  interpolated$interpolated <- FALSE
  
  i <- b <- 1
  while (i <= nrow(interpolated)) {
    if (!all(interpolated[i,1:2] == bedgraph[b,1:2])) {
      interpolated[i,5] <- TRUE
    } else {
      b <- b + 1
    }
    i <- i + 1
  }
  
  return(interpolated)
}

convertSlideToNonoverlapping <- function(bed) {
  return(interpolateBedScores(bed, interpolate=FALSE))
}

interpolateBedScores <- function(bed, mark.interpolated=F, interpolate=TRUE) {
  ## Expects a 4-item data frame:  chr, start, end, score
  ## Must have a consistent step size
  colnames(bed) <- c("chr","start","end","score")
  
  step <- as.numeric(names(which.max(with(bed, table(start[-1] - start[-nrow(bed)])))))
  window <- with(bed[1,,drop=F], end-start)
  
  ## Process each chromosome separately
  processOneChr <- function(curr) {
    all.starts <- data.frame(start=seq(curr$start[1], curr$start[nrow(curr)], step))
    curr <- merge(all.starts, curr, all.x=T)[,c(2,1,3,4)]
    curr$end <- curr$start + step
    curr$chr <- unique(curr$chr)[1]

    if (interpolate) {
      if (mark.interpolated) curr$interpolated <- is.na(curr$score)
      
      ## Linear interpolation
      fit <- approx(curr$score, n=length(curr$score))$y
      
      ## Nonlinear interpolation
      ##curr.ts <- ts(curr$score, f=4)
      ##fit <- ts(rowSums(tsSmooth(StructTS(curr.ts))[,-2]))
      curr$score <- fit
    }
    
    return(curr)
  }

  result <- tapply(1:nrow(bed), bed$chr, function(indices) {
    processOneChr(bed[indices,])
  })

  ## Remove items from the list that result from chr being a factor with
  ## no elements
  result <- result[which(!unlist(lapply(result,is.null)))]
    
  final <- result[[1]]
  final <- final[1:nrow(final),]

  if (length(result) > 1) {
    for (i in 2:length(result)) final <- rbind(final, result[[i]])
  }

  if (mark.interpolated) {
    final$score.original <- final$score
    final$score.original[final$interpolated] <- NA
    final$score.interpolated <- final$score
    final$score.interpolated[!final$interpolated] <- NA
  }
  
  return(final)
}


buildMaskRle <- function(fullBed, maskedBed) {
  ## full bed = slide and count file with every possible window
  ## masked bed = slide and count / bed file that contains only windows that pass a given mask (generated with ScriptureV2)
  ## Returns an Rle where 1 = include this line in full bed (not masked), 0 = remove line from full bed (masked)
  
  fullCon <- file(fullBed, open='r')
  maskedCon <- file(maskedBed, open='r')

  include.rle <- Rle(numeric())
  buffered <- numeric()
  
  maskedEmpty <- FALSE
  maskedWindow <- NULL
  
  advanceMasked <- function() {
    if (!maskedEmpty & length(maskedLine <- readLines(maskedCon, n=1, warn=F)) > 0) {
      maskedWindow <<- strsplit(maskedLine,'\t')[[1]][1:3]
    } else {
      maskedEmpty <<- TRUE
      maskedWindow <<- NULL
    }
  }

  advanceMasked()
  
  while (length(fullLine <- readLines(fullCon, n=1, warn=F)) > 0) {
    fullWindow <- strsplit(fullLine,'\t')[[1]][1:3]
    if (!is.null(maskedWindow) & all(fullWindow == maskedWindow)) {
      buffered <- c(buffered,1)
      advanceMasked()
    } else {
      buffered <- c(buffered,0)
    }

    if (length(buffered) == 100000) {
      include.rle <- c(include.rle, Rle(buffered))
      buffered <- numeric()
    }
  }
  include.rle <- c(include.rle, Rle(buffered))

  close(fullCon)
  close(maskedCon)
  
  return(include.rle)
}



################################################################################
makeTransparent<-function(someColor, alpha=100) {
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
  

selectEveryNth <- function(x, n, start=1) {
  return(x[seq(start, length(x), n)])
}


################################################################################
## Aggregate gene analysis
## March 17, 2013

plotGeneTrace <- function(data, tp=NULL, data.shuffled=NULL, ylim=c(0,25), column="ratio", legend=F, axis=T) {
  if (!is.null(tp)) curr <- subset(data, time == tp) else curr <- data

  if (!is.null(data.shuffled)) shuffled <- subset(data.shuffled, time == tp & !mask)
  active.genes <- subset(curr, expressed & !mask)
  inactive.genes <- subset(curr, !expressed & !mask)

  shapeTrace <- function(x) {
    trace <- with(x, tapply(get(column), coordinate, mean))
    stderr <- with(x, tapply(get(column), coordinate, function(y) sd(y)/sqrt(length(y))))
    result <- data.frame(coordinate=as.numeric(names(trace)), enrichment=trace, stderr=stderr)
    return(result)
  }
  active.trace <- shapeTrace(active.genes)
  inactive.trace <- shapeTrace(inactive.genes)
  if (!is.null(data.shuffled)) shuffled.trace <- shapeTrace(shuffled)

  with(active.trace, {
	plot(coordinate, enrichment, type='l', ylim=ylim, ylab='', xlab='', xaxt='n', yaxt='n', frame=F)
	polygon(c(coordinate, rev(coordinate)), c(enrichment + stderr*2, rev(enrichment - stderr*2)), col="#FF000044", border=NA)
	lines(coordinate, enrichment, type='l', lwd=1, col='red')
  })

  if (!is.null(data.shuffled)) {
         with(shuffled.trace, {
	  polygon(c(coordinate, rev(coordinate)), c(enrichment + stderr*2, rev(enrichment - stderr*2)), col="#00000044", border=NA)
	  lines(coordinate, enrichment, type='l', col='black', lwd=1)
	 })
  }

  with(inactive.trace, {
	polygon(c(coordinate, rev(coordinate)), c(enrichment + stderr*2, rev(enrichment - stderr*2)), col="#0000FF44", border=NA)
	lines(coordinate, enrichment, type='l', col='blue', lwd=1)
  })

  abline(v=0, lty=2)
  abline(v=40000,lty=2)    
  axis(2, tck=0.02, las=1)
  #axis(3, at=c(-100000,-50000,0,40000,90000,140000), labels=c("-100kb","-50kb","TSS","TES","+50kb","+100kb"), tick=F) #tck=0.02, tick=F)
  if (axis) axis(3, at=c(-100000,0,40000,140000), labels=c("-100kb","TSS","TES","+100kb"), tck=0.02)
  mtext(tp, side=2, line=1, cex=0.7)

  if (legend) {
    legend.names <- c("Active genes","Inactive genes","Random regions","95% CI")
    legend.colors <- c("red","blue","black","light gray")
    
    if (is.null(data.shuffled)) {
      legend.names <- legend.names[-3]
      legend.colors <- legend.colors[-3]
    }
    legend("topright",legend.names,col=legend.colors,lwd=2,bty='n')
  }
}


loadAggregateRegions <- function(file, nextgen=FALSE, score.column=8) {
  ## This function defaults for PlotAggregateRegions from Scripture version 1
  ## Set nextgen to TRUE for ScriptureV2
  
  x <- read.delim(file, header=F) #, stringsAsFactors=F)

  if (!nextgen) {
    colnames(x)[1:8] <- c("gene","section","coordinate","chr","start","end","strand","ratio")
  } else {
    x <- x[,c(1:6,score.column,18,22)]
    colnames(x) <- c("region","section","coordinate","chr","start","end","score","target.count","control.count")
  }
  return(x)
}


setupAggregateRegions <- function(file, name, ES.expression=NULL, input.cutoff=0, chr=NULL, mask=NULL) {
  x <- loadAggregateRegions(file)
  x <- x[,1:15]
  x$time <- name

  if (!is.null(ES.expression)) {
    x <- subset(x, x$gene %in% ES.expression$refseq)
    x$expressed <- x$gene %in% subset(ES.expression, expressed)$refseq
    x <- x[order(x$gene),]
  }

  if (!is.null(chr)) {
    x <- x[x$chr == chr,]
  }
  
  if (!is.null(mask)) {
    x$mask <- x$start <= mask[1] & x$end >= mask[2]
  } else {
    x$mask <- FALSE
  }

  x <- subset(x, V11 >= input.cutoff)
  #x$normRatio <- unlist(with(x, tapply(ratio, gene, function(ratios) ratios / mean(ratios))))
  return(x)
}




## Plot results from broad.pda.rap.PlotAggregateRegions
plotAggregateRegionsHelper <- function(aggregateStatistics, add=T, ...) {
  if ("col" %in% names(list(...))) {
    col <- list(...)[["col"]]
  } else {
    col <- "black"
  }
  
  with(aggregateStatistics, {
    if (!add) plot(coordinate, score, type='l', ...)
    else lines(coordinate, score, type='l', ...)
    polygon(c(coordinate, rev(coordinate)), c(score + stderr*2, rev(score - stderr*2)), col=makeTransparent(col), border=NA)
    lines(coordinate, score, type='l', ...)
  })
}


aggregateRegionsToStatistics <- function(data, column, fn = mean) {
  trace <- tapply(data[,column], data$coordinate, fn)
  stderr <- tapply(data[,column], data$coordinate, function(x) sd(x)/sqrt(length(x)))
  result <- data.frame(coordinate=as.numeric(names(trace)), score=trace, stderr=stderr)
  return(result)
}


plotAggregateRegionsHeatmap <- function(data, score.column, separator.columns=NULL, sort.coordinates=NULL, sort=TRUE, region.order=NULL, label.rows=FALSE, colors=c("white","orange","red","black"), ...) {
  ## Returns the order of the original data rows 
  require(gplots)

  regions <- unique(data$region)
  coordinates <- sort(unique(data$coordinate))
  data$coordinate <- as.numeric(as.matrix(factor(data$coordinate, levels=coordinates, labels=1:length(coordinates))))

  data.mat <- do.call(rbind, tapply(1:nrow(data), data$region, function(i) {
    scores <- rep(NA, length(coordinates))
    scores[data$coordinate[i]] <- data[i,score.column]
    scores
  }))
  colnames(data.mat) <- coordinates

  n <- apply(data.mat, 1, function(row) sum(is.na(row)))
  data.mat <- data.mat[n < length(coordinates)*0.4,]

  o <- 1:nrow(data.mat)
  if (!is.null(region.order)) {
    o <- region.order
    data.mat <- data.mat[o,]
  } else if (!is.null(sort.coordinates) & sort) {
    sort.value <- apply(data.mat, 1, function(row) mean(row[sort.coordinates], na.rm=T))
    o <- order(sort.value, decreasing=T)
    data.mat <- data.mat[o,]
  } else if (sort) {
    ## What is this for?
    require(vegan)
    d <- vegdist(data.mat, metric="gower", na.rm=T)
    h <- hclust(d)
    d <- as.dendrogram(h)
  }

  if (!is.null(separator.columns)) {
    for (i in 1:length(separator.columns)) {
      data.mat <- insertColumn(data.mat, rep(NA, nrow(data.mat)), separator.columns[i] + i - 1)
    }
  }
  
  heatmap.2(data.mat, Rowv=NA, Colv=NA, col = colorRampPalette(colors, space="rgb")(50), scale="none", labRow=if(label.rows) NULL else NA, labCol=NA, dendrogram="none", trace="none", density.info="none", ...)

  if (!is.null(region.order)) {
    return(region.order)
  } else if (!is.null(sort.coordinates) & sort) {
    return(o)
  } else if (sort) {
    return(h$order)
  } else {
    return(o)
  }
}


insertColumn <- function(mat, new.data, index) {
  ## Insert data into new column at given index
  if (index == 1) {
    begin <- c()
  } else {
    begin <- 1:(index-1)
  }

  if (index == ncol(mat) + 1) {
    end <- c()
  } else {
    end <- index:ncol(mat)
  }

  y <- mat[,c(begin,1,end)]
  y[,index] <- new.data
  return(y)
}


################################################################################
## Manhattan plots

manhattan <- function(...) {
  ## Looks nice for width=6, height=3
  par(cex=0.7, cex.lab=0.9, cex.axis=0.9, cex.main=0.5, mgp=c(3,0.5,0), oma=rep(0,4), mar=c(4.5,4.5,1,1), tck=-0.02)
  manhattan.helper(...) 
}

manhattan.helper <- function(bed, column="score", colors=c("gray10", "gray50"), ymax="max", cutoffline=10, annotate=NULL, downsample=NULL, ylab=expression(-log[10](italic(p))), ...) {  
  bed <- bed[with(bed, order(chr, start, end)),]

  if (!is.null(downsample)) {
    ## Downsample points below the given cutoff
    below <- which(bed[,column] < downsample)
    bed <- bed[-1*sample(below, floor(length(below)*0.95)),]
  }

  bed[is.infinite(bed[,column]) & sign(bed[,column]) == 1,column] <- 350
  
  bed$mid <- with(bed, floor((start + end) / 2))
  nchr <- length(unique(bed$chr))
  bed$pos = NA
  ticks=NULL
  lastbase=0
  colors <- rep(colors,nchr)[1:nchr]
  if (ymax=="max") ymax<-ceiling(max(bed[,column], na.rm=T))
    
  if (nchr == 1) {
    bed$pos = bed$mid
    ticks = floor(length(bed$pos))/2+1
  } else {
    for (i in 1:nchr) {
      chr <- unique(bed$chr)[i]
      which.chr <- which(bed$chr == chr)
      if (i==1) {
        bed[which.chr,]$pos = bed[which.chr, ]$mid
      } else {
        bed[which.chr,]$pos <- lastbase + bed[which.chr,]$mid
      }
      ticks <- c(ticks, bed[which.chr,]$pos[floor(length(bed[which.chr, ]$pos)/2)+1])
      lastbase <- lastbase + bed$mid[tail(which.chr,1)]
    }
  }
    
  if (nchr==1) {
    with(bed, plot(pos, get(column), ylim=c(0,ymax), pch=19, cex=0.5, ylab=ylab, ...)) #xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  } else {
    with(bed, plot(pos, get(column), ylim=c(0,ymax), xaxt='n', type='n', xlab="Chromosome", ylab=ylab, ...)) #xlab="Chromosome", xaxt="n", type="n", ...))
    axis(1, at=ticks, lab=unique(bed$chr), las=2, ...)
    icol=1
    for (i in unique(bed$chr)) {
      with(subset(bed, chr == i), lines(pos, get(column), col=colors[icol], pch=19, cex=0.5, ...))
      icol=icol+1
    }
  }
  
  #if (!is.null(annotate)) {
  #  d.annotate=bed[which(d$SNP %in% annotate), ]
  #  with(d.annotate, points(pos, logp, col="green3", ...)) 
  #}
  
  if (cutoffline) abline(h=cutoffline, col="red")
}







################################################################################
## GViz utils

separateTrack <- function(data, mask) {
  ## separate a track into two parts, filling in the remainder with NA
  result <- rbind(data, data)
  result[1,mask] <- NA
  result[2,!mask] <- NA
  return(result)
}

getRefseqTrack <- function(chr, from, to, ...) {
  refseq <- UcscTrack(genome="mm9", chromosome=as.character(chr), track="refGene", from=from, to=to, trackType = "GeneRegionTrack",
                    rstart = "exonStarts", rends = "exonEnds", gene = "name2", symbol = "name2", transcript = "name",
                    strand = "strand", name = "Genes", fill="blue", collapseTranscripts=TRUE, lwd=0, showId=TRUE, ...)
  return(refseq)
}

interpolated.colors <- c("#666666", "#444444") #c("#A7A9AC","#58595B")
score.names <- c("score.interpolated","score.original")
makeDataTrack <- function(bed, chr=NULL, interpolate=T, log2=F, center=F, smooth=NA, separate=F, ...) {

  if (!is.null(chr)) {
    bed <- bed[bed$chr == chr,]
  }

  if (log2) {
    bed$score <- log2(bed$score)
  }

  if (center) {
    bed$score <- scale(bed$score, center=T, scale=F)[,1]
  }

  if (!is.na(smooth)) {
    bed$score <- filter(bed$score, rep(1/smooth, smooth))
  }

  args <- list(...)
  if (is.null(args$type)) args$type <- 'h'
  if (is.null(args$lwd)) args$lwd <- 0.25

  if (!interpolate) {

    if (separate) {
      st <- separateTrack(bed$score, bed$score >= 0)
      rownames(st) <- c("positive","negative")
      args$data <- st
      if (is.null(args$col)) args$col <- c("red","blue")
      args$groups <- rownames(st)
    } else {
      args$data <- bed$score
    }
    
  } else {
    bed <- interpolateBedScores(bed[,c("chr","start","end","score")], mark.interpolated=T)
    args$data <- bed[,score.names]
    args$groups <- score.names
    if (is.null(args$col)) args$col <- interpolated.colors
  }

  args$start <- bed$start
  args$end <- bed$end
  args$chromosome <- as.character(bed$chr)

  return(do.call(DataTrack, args))
}

loadMm9Sizes <- function() {
  sizes <- read.delim("/seq/lincRNA/data/mm9/sizes", header=F, stringsAsFactors=F)
  x <- sizes$V2
  names(x) <- sizes$V1
  return(x)
}

getAllCoverageFromBam <- function(bamFile, sizes=loadMm9Sizes(), ...) {
  require(rtracklayer)
  require(Rsamtools)
  result <- RleList()
  for (chr in names(sizes)) {
    result[[chr]] <- getCoverageFromBam(chr, 0, sizes[chr], ...)
    #param <- ScanBamParam(what = c("rname", "pos", "qwidth", "strand"),
    #                      which = GRanges(chr, IRanges(0,sizes[chr])),
    #                      flag = scanBamFlag(isUnmappedQuery = FALSE))
    #result[[chr]] <- getCoverageFromBamHelper(bamFile, param, width=sizes[chr], ...)
  }
  return(result)
}

getCoverageFromBamHelper <- function(bamFile, param, midpoints=FALSE, width=NULL, collapseTo5p=FALSE, extend=NA, ...) {
  require(Rsamtools)
  require(Gviz)
  
  x <- scanBam(bamFile, ..., param = param)[[1]]
  if (length(x$pos) == 0) {
    cvg <- Rle(0,width)
  } else if (midpoints) {
    cvg <- coverage(IRanges(x[["pos"]] + x[["qwidth"]]/2, width = 1), width=width)
  } else if (collapseTo5p) {
    cvg <- coverage(IRanges(x[["pos"]][x$strand == "+"], width = 1), width=width) +
      coverage(IRanges((x[["pos"]] + x[["qwidth"]])[x$strand == "-"], width=1), width=width)
  } else if (!is.na(extend)) {
    cvg <- coverage(IRanges(x[["pos"]][x$strand == "+"], width=extend), width=width) + 
      coverage(IRanges((x[["pos"]] + x[["qwidth"]])[x$strand == "-"], width=extend), width=width)
  } else {
    cvg <- coverage(IRanges(x[["pos"]], width = x[["qwidth"]]), width=width)
  }
  return(cvg)
}

getCoverageFromBam <- function(chr, start, end, bamFile, strand=NA, matchFirstInPairStrand=NA, ...) {
  require(Rsamtools)
  if (!is.na(strand)) {
    if (is.na(matchFirstInPairStrand)) {
      flag <- scanBamFlag(isUnmappedQuery=FALSE,
                          isMinusStrand=(strand=='-'))
    } else if (matchFirstInPairStrand == "same") {
      flag <- scanBamFlag(isUnmappedQuery=FALSE,
                          isFirstMateRead=TRUE,
                          isMinusStrand=(strand=='-'))
    } else if (matchFirstInPairStrand == "opposite") {
      flag <- scanBamFlag(isUnmappedQuery=FALSE,
                          isFirstMateRead=TRUE,
                          isMinusStrand=(strand!='-'))
    } else {
      warning("Unexpected value for matchFirstInPairStrand")
    }
  } else {
    flag <- scanBamFlag(isUnmappedQuery=FALSE)
  }
  param <- ScanBamParam(what = c("pos", "qwidth", "strand"),
                        which = GRanges(chr, IRanges(start, end)),
                        flag = flag)
  return(getCoverageFromBamHelper(bamFile, param, width=end, ...))
}

getCoverageFromBED <- function(chr, start, end, bedFile, strand=NA, ...) {
  require(Gviz)
  require(GenomicRanges)
  bed <- readBed(bedFile)
  bed <- bed[bed$chr == chr,]
  if (!is.na(strand)) bed <- bed[bed$strand == strand,]
  cvg <- with(bed, coverage(IRanges(start, width=end-start)))
  return(cvg)
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

getCoverageTrackFromBam <- function(chr, start, end, bamFile, res=10, overlap=0, FUN=mean, percent=FALSE, midpoints=FALSE, collapseTo5p=FALSE, strand=NA, matchFirstInPairStrand=NA, ...) {
  ## If percent is true, then res is interpreted to be a percentage
  ## (e.g., res=10 will give 10 windows each containing 10% of start to end)
  
  cvg <- getCoverageFromBam(chr, start, end, bamFile, midpoints=midpoints, collapseTo5p=collapseTo5p, strand=strand, matchFirstInPairStrand=matchFirstInPairStrand)
  return(getTrackFromCoverage(cvg, chr, start, end, res=res, overlap=overlap, FUN=FUN, percent=percent, ...))
}


getCoverageTrackFromBED <- function(chr, start, end, bedFile, res=10, overlap=0, FUN=mean, percent=FALSE, midpoints=FALSE, collapseTo5p=FALSE, strand=NA, ...) {
  cvg <- getCoverageFromBED(chr, start, end, bamFile, midpoints=midpoints, collapseTo5p=collapseTo5p, strand=strand)
  return(getTrackFromCoverage(cvg, chr, start, end, res=res, overlap=overlap, FUN=FUN, percent=percent, ...))
}


getTrackFromCoverage <- function(cvg, chr, start, end, res=10, overlap=0, FUN=mean, percent=FALSE, ...) {
  if (percent) {
    intervals <- seq(start, end, length.out=(ceiling(100 / res) + 1))
    windows <- IRanges(start=intervals[-length(intervals)], end=intervals[-1])
  } else {
    windows <- IRanges(start=seq(start, end-res, res-overlap), end=NA, width=res)
  }

  tryCatch({
    scores <<- IRanges:::aggregate(cvg, windows, FUN=FUN)
  }, error = function(e) {
    scores <<- aggregate(cvg, windows, FUN=FUN) ## This works in R-3.1
  })

  args <- list(...)
  if (is.null(args$type)) args$type <- 'h'
  if (is.null(args$lwd)) args$lwd <- 0.25
  args$chromosome <- chr
  args$range <- windows
  args$data <- scores
  return(do.call(DataTrack, args))  
}


getEnrichmentTrackFromBams <- function(chr, start, end, bamFileNum, bamFileDenom, res=10, FUN=mean, offset=0.001, relative.offset=TRUE, num.total=NA, den.total=NA, collapseTo5p=FALSE, ...) {
  ## offset = raw reads to add to each count (relative.offset = FALSE) or fragments per million to add to each count (relative.offset = TRUE)
  
  if (is.na(num.total)) num.total <- countBamMapped(bamFileNum)
  if (is.na(den.total)) den.total <- countBamMapped(bamFileDenom)

  numerator <- getCoverageTrackFromBam(chr, start, end, bamFileNum, res=res, FUN=FUN, collapseTo5p=collapseTo5p, ...)
  denominator <- getCoverageTrackFromBam(chr, start, end, bamFileDenom, res=res, FUN=FUN, collapseTo5p=collapseTo5p, ...)

  args <- list(...)
    if (is.null(args$type)) args$type <- 'h'
  if (is.null(args$lwd)) args$lwd <- 0.25
  args$chromosome <- chr
  args$range <- numerator@range
  data <- matrix(NA, nrow=3, ncol=length(numerator@data), dimnames=list(rownames=c("ratio","count.target","count.control")))

  if (relative.offset) {
    num.offset <- num.total / 1000000 * offset
    den.offset <- den.total / 1000000 * offset
  } else {
    num.offset <- offset
    den.offset <- offset
  }
  
  data[1,] <- (numerator@data + num.offset) / (denominator@data + den.offset) * (den.total / num.total)
  data[2,] <- numerator@data
  data[3,] <- denominator@data
  args$data <- data
  return(do.call(DataTrack, args))  
}


countBamMapped <- function(bamFile) {
  ## Counts number of reads in a BAM file WITHOUT iterating.  Requires that the BAM is indexed
  y <- do.call(rbind, strsplit(system(paste("samtools idxstats",bamFile), intern=TRUE),"\t"))
  return(sum(as.numeric(y[-nrow(y),3])))
}

countBamTotal <- function(bamFile) {
  ## Counts number of reads in a BAM file WITHOUT iterating.  Requires that the BAM is indexed
  y <- do.call(rbind, strsplit(system(paste("samtools idxstats",bamFile), intern=TRUE),"\t"))
  return(sum(as.numeric(y[,4])))
}


findOverlapsOppositeStrand <- function(granges1, granges2, ...) {
  overlaps <- findOverlaps(granges1, granges2, ignore.strand=TRUE, ...)
  to.keep <- granges1[queryHits(overlaps)]@strand != granges2[subjectHits(overlaps)]@strand
  return(overlaps[to.keep])
}
  

myUcscTrack <- function (track, table = NULL, trackType = c("AnnotationTrack", 
    "GeneRegionTrack", "DataTrack", "GenomeAxisTrack"), genome, 
    chromosome, name = NULL, from, to, ...) 
{
    trackType <- match.arg(trackType)
    if (missing(genome) || !IRanges:::isSingleString(genome)) 
        stop("Need to specify genome for creating a UcscTrack")
    if (missing(chromosome)) 
        stop("Need to specify chromosome for creating a UcscTrack")
    chromosome <- .chrName(chromosome)[1]
    sessionInfo <- Gviz:::.cacheTracks(genome = genome, chromosome = chromosome, 
        track = track, env = Gviz:::.ucscCache)
    if (missing(from)) 
        from <- 1
    if (missing(to)) 
        to <- sessionInfo$chrInfo[chromosome]
    gr <- GRanges(ranges = IRanges(start = from, end = to), seqnames = chromosome)
    suppressWarnings(genome(gr) <- unname(genome))[1]
    query <- rtracklayer:::ucscTableQuery(sessionInfo$session, sessionInfo$track, 
        gr)
    if (!is.null(table)) {
        table <- match.arg(table, sessionInfo$availTables)
        rtracklayer:::tableName(query) <- table
    }
    if (is.null(name)) 
        name <- if (is.null(table)) 
            track
        else paste(sessionInfo$track, table)
    tableDat <- if (trackType == "DataTrack") {
        tmp <- try(rtracklayer:::track(query), silent = TRUE)
        if (is(tmp, "try-error")) {
            warning(tmp)
            data.frame()
        }
        else as.data.frame(tmp)
    }
    else {
        tmp <- try(rtracklayer:::getTable(query), silent = TRUE)
        if (is(tmp, "try-error")) {
            warning(tmp)
            data.frame()
        }
        else tmp
    }
    if (is(tmp, "try-error") && nrow(tableDat) == 0) 
        stop("Error fetching data from UCSC")
    args <- lapply(list(...), function(x) if (is.character(x) && 
        length(x) == 1) 
        if (!x %in% colnames(tableDat)) 
            x
        else tableDat[, x]
    else x)
    if (trackType == "GeneRegionTrack") {
        args$start <- from
        args$end <- to
    }
    args <- lapply(args, function(x) if (!length(x)) 
        NULL
    else x)
    trackObject <- do.call(trackType, args = c(list(chromosome = chromosome, 
        genome = genome, name = name), args))
    return(trackObject)
}



myPlotTracks <- function(tracks, from, to, ...) {
  p <- plotTracks(tracks, from=from, to=to, cex=0.5, cex.axis=0.5, cex.title=0.4, col.axis='black', background.title='white', col.title='black')
  return(p)
}


Rle2DataBits <- function(x) {
  sel <- suppressWarnings(runValue(x) != 0)
  list(start=start(x)[sel], end=end(x)[sel], data=runValue(x)[sel])
}


grAggregateScores <- function(gr, windows, FUN=mean, col="score", default.val=0) {
  ## Slide over genomic ranges and apply FUN to the specified column (containing numeric score data)
  ## Return genomic ranges windows with an added mcol called $FUN.$col
  ## Assumes all elements in gr are of the same size
  ## default.val = value to fill in when there are no overlapping elements in gr for a window

  grWidth <- with(gr, end(ranges(gr))-start(ranges(gr)))[1]
  overlaps <- as.data.frame(findOverlaps(gr, windows, ignore.strand=TRUE, minoverlap=grWidth))
  scores <- tapply(mcols(gr)[[col]][overlaps$queryHits], overlaps$subjectHits, FUN)
  result <- rep(default.val, length(windows))
  result[as.numeric(names(scores))] <- scores
  mcols(windows)[[paste0(col,".",as.character(substitute(FUN)))]] <- result
  return(windows)
}


strandedBamImport <- function (file, selection) { 
  ## Import and split reads by first-in-pair strand
  ## This isn't doing what I want it to do ...
  require(Rsamtools)
  if (!file.exists(paste(file, "bai", sep = "."))) stop("Unable to find index for BAM file '", file, "'. You can build an index using the following command:\n\t", "library(Rsamtools)\n\tindexBam(\"", file, "\")")
  sinfo <- scanBamHeader(file)[[1]]
  res <- if (!as.character(seqnames(selection)[1]) %in% names(sinfo$targets)) { 
    mcols(selection) <- DataFrame(score = 0)
    selection 
  } else { 
    param <- ScanBamParam(what = c("pos", "qwidth","strand","flag"), which = selection, flag = scanBamFlag(isUnmappedQuery = FALSE))
    x <- scanBam(file, param = param)[[1]]
    
    ## This splits by read strand:
    #gr <- GRanges(strand=x[["strand"]], ranges=IRanges(x[["pos"]], width = x[["qwidth"]]), seqnames=seqnames(selection)[1])
    ## This splits by first in pair strand
    gr <- GRanges(strand=c("-","+")[as.numeric(bamFlagTest(x[["flag"]], "isFirstMateRead") != (x[["strand"]]=="+"))+1], 
                  ranges=IRanges(x[["pos"]], width = x[["qwidth"]]), 
                  seqnames=seqnames(selection)[1])

    grs <- split(gr, strand(gr))
    cov <- lapply(grs[c("+", "-")], function(y) coverage(ranges(y), width=end(selection)))
    pos <- sort(unique(unlist(lapply(cov, function(y) c(start(y), end(y))))))
    if (length(pos)==0) { 
      mcols(selection) <- DataFrame(plus=0, minus=0)
      selection 
    } else { 
      GRanges(seqnames = seqnames(selection)[1], ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)), 
        plus=as.numeric(cov[["+"]][head(pos, -1)]), 
        minus=-as.numeric(cov[["-"]][head(pos, -1)])) 
    }
  } 
  return(res)
} 

makeDirectionalCTCFTrack <- function(bamFile, broadPeak.bed, motif.bed, chr, start, end, col=c("pink","gray","purple"), ...) {
  broadPeak <- readBed(broadPeak.bed)
  motifs <- readBed(motif.bed); colnames(motifs)[4:6] <- c("strand","score","p")

  ## Assign each broadPeak to one or both strands
  ixn <- findOverlaps(GRangesFromBed(motifs[,1:3]), GRangesFromBed(broadPeak[,1:3]))
  res <- do.call(rbind, tapply(queryHits(ixn), subjectHits(ixn), function(motif.i) {
    with(motifs[motif.i,], tapply(score, strand, sum))
  }))

  pos.peaks <- GRangesFromBed(broadPeak[as.numeric(rownames(res)[!is.na(res[,1]) & (res[,1] > 20 | is.na(res[,2]) | res[,1] > res[,2]*2)]),1:3])
  neg.peaks <- GRangesFromBed(broadPeak[as.numeric(rownames(res)[!is.na(res[,2]) & (res[,2] > 20 | is.na(res[,1]) | res[,2] > res[,1]*2)]),1:3])
  other.peaks <- GRangesFromBed(broadPeak[-as.numeric(rownames(res)),1:3])

  selection <- GRanges(seqnames=chr, ranges=IRanges(start,end))
  require(Rsamtools)
  if (!file.exists(paste(bamFile, "bai", sep = "."))) stop("Unable to find index for BAM file '", file, "'. You can build an index using the following command:\n\t", "library(Rsamtools)\n\tindexBam(\"", file, "\")")

  sinfo <- scanBamHeader(bamFile)[[1]]
  res <- if (!as.character(seqnames(selection)[1]) %in% names(sinfo$targets)) { 
    mcols(selection) <- DataFrame(score = 0)
    selection 
  } else { 
    param <- ScanBamParam(what = c("pos", "qwidth"), which = selection, flag = scanBamFlag(isUnmappedQuery = FALSE))
    x <- scanBam(bamFile, param = param)[[1]]

    gr <- GRanges(ranges=IRanges(x[["pos"]], width = x[["qwidth"]]), 
                  seqnames=seqnames(selection)[1])

    grs <- list(`+`=subsetByOverlaps(gr, pos.peaks),
                `-`=subsetByOverlaps(gr, neg.peaks),
                other=subsetByOverlaps(gr, other.peaks))
    cov <- lapply(grs[c("+", "-", "other")], function(y) coverage(ranges(y), width=end(selection)))
    pos <- sort(unique(unlist(lapply(cov, function(y) c(start(y), end(y))))))
    if (length(pos)==0) { 
      mcols(selection) <- DataFrame(plus=0, minus=0)
      selection 
    } else {
      GRanges(seqnames = seqnames(selection)[1], ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)), 
        plus=as.numeric(cov[["+"]][head(pos, -1)]), 
        minus=-as.numeric(cov[["-"]][head(pos, -1)]),
        other=as.numeric(cov[["other"]][head(pos, -1)])) 
    }
  }
  return(DataTrack(start=start(res), end=end(res), data=t(as.matrix(mcols(res))), chromosome=chr, type="h", col=col, groups=c("plus","minus","other"), ...))
}

#####################################################################
## Calculate and plot Z-scores
getLocalZscore <- function(row, bed, distance=5000000) {
  scores <- getBedIntersection(as.character(row$chr), row$start, row$end, bed, distance)$score
  return(c((row$score - mean(scores, na.rm=T)) / sd(scores, na.rm=T), (row$score / mean(scores, na.rm=T))))
}

findInitiationSites <- function(bed, smooth=1, distance=5000000, z.cutoff=2, log2=F) {
  if (length(unique(bed$chr)) > 1) stop("Takes forever with multiple chromosomes")
  bed$score <- filter(bed$score, rep(1/10,10))
  if (log2) bed$score <- log2(bed$score)
  z.fold <- sapply(1:nrow(bed), function(row) getLocalZscore(bed[row,,drop=F], bed, distance))
  z <- z.fold[1,]
  fold <- z.fold[2,]
  z.track <- with(bed, DataTrack(data=z, start=start, end=end, chromosome=as.character(chr),
                                 name="Local Z-score", type='l', col='black', lwd=0.75, size=3))
  if (z.cutoff > 0) z.initiation <- bed[!is.na(z) & z >= z.cutoff, 1:3]
  else z.initiation <- bed[!is.na(z) & z <= z.cutoff, 1:3]
  z.sites.track <- with(z.initiation, AnnotationTrack(start=start, end=end, chromosome=as.character(chr),
                                                    name="Z Sites", fill="black", stacking="dense", lwd=0))

  bed$fold <- fold
  bed <- merge(bed, z.initiation)
  mean.fold <- with(bed, mean(fold, na.rm=T))
  
  return(list(z=z, track=z.track, sites=z.initiation, sites.track=z.sites.track, fold=mean.fold))
}

getBedCorrelation <- function(bed1, bed2, col="score", col1=col, col2=col, ...) {
  ## Try to merge bed files based on chr, start, end, and calculate the correlation of their score columns
  x <- merge(bed1[,c("chr","start","end",col1)], bed2[,c("chr","start","end",col2)], by=c("chr","start","end"))

  if (col1 == col2) {
    col1 <- paste(col1,".x",sep='')
    col2 <- paste(col2,".y",sep='')
  }
  
  cor(x[,col1], x[,col2], use="pairwise.complete.obs", ...)
}

smoothAndPrune <- function(bed, smooth=10) {
  bed$score <- filter(bed$score, rep(1/smooth, smooth))
  bed <- bed[selectEveryNth(1:nrow(bed),smooth),]
  return(bed)
}


#####################
## HiC Utilities

anchorTrackFromHicMatrix <- function(matrix, chr, bin.size, anchor.start, anchor.end, fun=sum, limit=NULL) {
  ## Takes a HiC matrix (e.g., from Erez's HiC Viewer), and outputs a bedgraph file for the given anchor
  bin.start <- ceiling(anchor.start / bin.size)
  bin.end <- ceiling(anchor.end / bin.size)

  if (bin.start == bin.end) {
    scores <- matrix[bin.start,]
  } else {
    scores <- apply(matrix[,bin.start:bin.end], 1, fun)
  }

  bins <- seq(0,nrow(matrix)*bin.size,bin.size)
  result <- data.frame(chr=chr, start=bins, end=(bin.size + bins), score=scores)
  return(result)
}

anchorTrackFromObservedHicMatrix <- function(matrix, chr, bin.size, anchor.start, anchor.end, normalization=NULL, fun=sum, limit=NULL) {
  ## HiC Viewer observed matrices have col 1: bin 1, col 2: bin 2, col 3: score / count
  ## and it is not symmetric - gives only one triangle half of the matrix

  ## Normalization table info:
  ## The three *norm files are normalization vectors that can be used to transform the raw contact matrices M into normalized matrices M*. (See the glossary below and section II.b of the Extended Experimental Procedures of Rao, Huntley, et al, Cell, 2014 for more information about the different types of normalizations. All analyses and results presented in the main text of Rao, Huntley, et al, Cell, 2014 were performed using KR normalized contact matrices.) Each file is ordered such that the first line of the normalization vector file is the norm factor for the first row/column of the corresponding raw contact matrix, the second line is the factor for the second row/column of the contact matrix, and so on. To normalize, an entry M_i,j in a *RAWobserved file, divide the entry by the corresponding norm factors for i and j. 

  ## For example, here is a line from the 5kb chr1 MAPQGE30 raw observed contact matrix (K562/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.RAWobserved):
  ## 40000000  40100000  13.0

  ## To normalize this entry using the KR normalization vector, one would divide 59.0 by the 8001st line ((40000000/5000)+1=8001) and the 8021st line ((40100000/5000)+1=8021) of K562/5kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_5kb.KRnorm. The 8001st line of the KR norm file is 1.391205327929388;The 8021st line of the KR norm file is 1.764507526190848. So the corresponding KR normalized entry for the entry above is 13.0/(1.391205327929388*1.764507526190848) or 5.2957637802416944. If the KR normalization vector file is empty or all NaNs, then the KR algorithm didn’t converge on that particular matrix (likely due to sparsity of the matrix).  

  bin.start <- floor(anchor.start / bin.size)*bin.size
  bin.end <- floor(anchor.end / bin.size)*bin.size

  y <- subset(matrix, V1 >= bin.start & V1 <= bin.end)
  if (!is.null(limit)) y <- subset(y, V2 >= bin.start - limit & V2 <= bin.end + limit)
  if (!is.null(normalization)) {
    y$V3 <- apply(y, 1, function(row) row[3]/(normalization[(row[1]/bin.size)+1]*normalization[(row[2]/bin.size)+1]))
  }

  scores <- with(y, tapply(V3, V2, sum))
  result <- data.frame(chr=chr, start=as.numeric(names(scores)), end=as.numeric(names(scores)) + bin.size, score=scores)

  y <- subset(matrix, V2 >= bin.start & V2 <= bin.end)
  if (!is.null(limit)) y <- subset(y, V1 >= bin.start - limit & V1 <= bin.end + limit)
  if (!is.null(normalization)) {
    y$V3 <- apply(y, 1, function(row) row[3]/(normalization[(row[1]/bin.size)+1]*normalization[(row[2]/bin.size)+1]))
  }
  scores <- with(y, tapply(V3, V1, sum))
  result2 <- data.frame(chr=chr, start=as.numeric(names(scores)), end=as.numeric(names(scores)) + bin.size, score=scores)

  result <- rbind(result2, result)
  result <- result[order(result$start),]
  result <- subset(result, !duplicated(start)) # not sure if this is the right solution here ...  fix better?
  return(result)
}


getChromosomeSize <- function(chr, genome="mm9") {
  if (genome != "mm9") return -1;
  sizes <- read.delim("/seq/lincRNA/data/mm9/sizes", header=F)
  return(subset(sizes, V1 == chr)$V2)
}


#####################
## Methylation Utils

loadBisSNPCalls <- function(sample.name, basename) {
  ## Load CpG coverage and methylation data from bedgraph files output by BisSNP  
  methylated <- read.delim(paste(basename, ".cpg.filtered.sort.CG.bedgraph", sep=''), header=F)
  coverage <- read.delim(paste(basename, ".cpg.filtered.sort.CG.coverage.bedgraph", sep=''), header=F)
  methylated$V5 <- coverage$V4
  colnames(methylated) <- c("chr","start","end","pct","coverage")
  methylated$reads <- with(methylated, round(coverage*pct / 100))
  methylated$sample <- sample.name
  return(methylated)
}


loadMultipleBisSNPCalls <- function(sample.name, ...) {
  ## Load CpG coverage and methylation data from bedgraph files output by BisSNP  
  ## arguments to ... should be the basename of the BisSNP data
  basenames <- list(...)

  all.data <- loadBisSNPCalls(sample.name, basenames[[1]])
  for (basename in unlist(basenames)[-1]) {
    curr <- loadBisSNPCalls(sample.name, basename)
    all.data <- merge(all.data, curr, by=c("chr","start","end","sample"), all=T)
    all.data$coverage.x[is.na(all.data$coverage.x)] <- 0
    all.data$coverage.y[is.na(all.data$coverage.y)] <- 0
    all.data$reads.x[is.na(all.data$reads.x)] <- 0
    all.data$reads.y[is.na(all.data$reads.y)] <- 0
    all.data$coverage <- with(all.data, coverage.x + coverage.y)
    all.data$reads <- with(all.data, reads.x + reads.y)
    all.data <- all.data[,c("chr","start","end","coverage","reads","sample")]
  }
  all.data$pct <- with(all.data, reads / coverage * 100)
  all.data <- all.data[,c("chr","start","end","pct","coverage","reads","sample")]
  return(all.data)
}


getVectorRevComp <- function(seqs) {
  require(Biostrings)
  sapply(seqs, function(x) as.character(reverseComplement(DNAString(x))))
}



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  require(grid)
  
                                        # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
                                        # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
                                        # Make the panel
                                        # ncol: Number of columns of plots
                                        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
      }

  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
                                        # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

                                        # Make each plot, in the correct location
    for (i in 1:numPlots) {
                                        # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                          layout.pos.col = matchidx$col))
    }
  }
}


divideMatrixByCol <- function(mat, col.vals) {
  if (nrow(mat) != length(col.vals)) cat("Dimensions are incompatible\n")
  return(mat / col.vals)
}

divideMatrixByRow <- function(mat, row.vals) {
  if (ncol(mat) != length(row.vals)) cat("Dimensions are incompatible\n")
  return(t(t(mat) / row.vals))
}


divideRowsByValues <- function(mat, row.vals) {
  return(divideMatrixByCol(mat, row.vals))
}


divideColsByValues <- function(mat, col.vals) {
  return(divideMatrixByRow(mat, col.vals))
}



###################################################################################
## August 27, 2015
## TSS analysis

getBestTSS <- function(tss) {
  ## tss: BED file, where name field has format Gene_TSSX and score is the % explained
  ## output: subset of the original BED file, filtered to only the TSS per gene with the highest score

  tss$gene <- sapply(as.character(as.matrix(tss$name)), function(s) strsplit(s, "_TSS")[[1]][1])
  mask <- na.omit(unlist(with(tss, tapply(1:nrow(tss), gene, function(indices) indices[which.max(score[indices])]))))
  return(tss[mask,])
}

getSecondBestTSS <- function(tss, cutoff=0.5) {
  ## tss: BED file, where name field has format Gene_TSSX and score is the % explained
  ## cutoff:  get the second best TSS but only if it has >cutoff signal compared to the best TSS
  ## output: subset of the original BED file, filtered to only the single second best TSS per gene with the highest score, if one passes the cutoff
  tss$gene <- sapply(as.character(as.matrix(tss$name)), function(s) strsplit(s, "_TSS")[[1]][1])
  mask <- na.omit(unlist(with(tss, tapply(1:nrow(tss), gene, function(indices) {
    indices <- indices[order(score[indices], decreasing=T)]
    if (length(indices) > 1 & score[indices[2]] >= 0.5 * score[indices[1]]) {
      return(indices[2])
    } else {
      return(c())
    }
  }))))
  return(tss[mask,])
}

getBestAndSecondBestTSS <- function(tss, cutoff=0.5) {
  return(rbind(getBestTSS(tss), getSecondBestTSS(tss, cutoff)))
}




getGeneNameMap <- function(genome="mm9", write.file=FALSE) {
  if (genome == "mm9") {
    dataset='mmusculus_gene_ensembl'
    host="may2012.archive.ensembl.org"
    symbol.col="mgi_symbol"
    file="/seq/lincRNA/data/mm9/mm9.refseqToSymbol.bed"
  } else if (genome == "hg19") {
    dataset='hsapiens_gene_ensembl'
    host="may2012.archive.ensembl.org"
    symbol.col="hgnc_symbol"
    file="/seq/lincRNA/data/hg19/hg19.refseqToSymbol.bed"
  } else {
    stop(paste0("Genome ",genome," not recognized"))
  }

  if (file.exists(file) & !write.file) {
    genes <- read.delim(file, stringsAsFactors=F)
    return(genes)
  } else {
    library(biomaRt)
    attributes=c("ensembl_gene_id", "refseq_mrna", "refseq_ncrna", symbol.col, "chromosome_name", "start_position", "end_position", "strand", "transcript_start", "transcript_end")
    mart <- useMart('ENSEMBL_MART_ENSEMBL',dataset=dataset, host=host)
    all.genes <- getBM(attributes=attributes, mart=mart)
    genes <- with(all.genes, data.frame(chr = paste0("chr", chromosome_name), start=transcript_start, end=transcript_end,
                                        name=ensembl_gene_id,
                                        score=0,
                                        strand = c("+","-")[as.numeric(strand < 0) + 1],
                                        symbol = get(symbol.col),
                                        refseq = refseq_mrna))
    genes$refseq <- as.character(as.matrix(genes$refseq))
    genes$refseq[genes$refseq == ""] <- as.character(as.matrix(all.genes$refseq_ncrna[genes$refseq == ""]))
    if (write.file) write.table(genes, file=file, sep='\t', quote=F, col.names=T, row.names=F)
    return(genes)
  }
}


getSimpleGeneNameMap <- function(genome="mm9", write.file=FALSE) {
  ## December 10, 2015
  ## including both mm9 protein-coding and lincRNA genes, one isoform per gene, just genomic span
  ## Chooses between multiple refseq isoforms by picking the lowest number (usually the dominant isoform)

  ## Added more genomes April 22, 2016
  if (genome == "mm9") {
    file="/seq/lincRNA/data/mm9/mm9.ESCGeneSetWithLncRNA.txt"
    linc.file="/seq/lincRNA/Annotations/lincRNA/140807_lncRNA.bed"
  } else if (genome == "hg19") {
    file="/seq/lincRNA/data/hg19/hg19.SimpleGeneNameMap.txt"
    linc.file=NULL
  } else {
    stop(paste0("Genome ",genome," not recognized"))
  }

  if (file.exists(file) & !write.file) {
    genes <- read.delim(file, stringsAsFactors=F)
    return(genes)
  } else {
    #linc.tss <- readBed("/seq/lincRNA/Jesse/CRISPR_Screen/150304_CRISPRi_Pool_Design/TSS/140807_lncRNA.TSS_100bp.bed")

    gene.map <- getGeneNameMap(genome=genome)
    gene.map$NM <- grepl("NM", gene.map$refseq)
    gene.map$refseqnum <- sapply(as.character(as.matrix(gene.map$refseq)), function(x) strsplit(x, "_")[[1]][2])

    gene.map <- subset(gene.map, symbol != "")
    to.use <- with(gene.map, tapply(1:nrow(gene.map), symbol, function(i) {
      ## Prefer things with refseq entries, or 
      if (any(NM[i])) {
        i <- i[which(NM[i])]
      } else if (any(refseq[i] != "")) {
        i <- i[which(refseq[i] != "")]
      }
      result <- i[which.min(refseqnum[i])[1]]
      if (length(result) == 0 | is.na(result)) 
        result <- i[1]
      return(result)
      }))
    unique.symbol <- gene.map[to.use,]
    unique.symbol$name <- unique.symbol$symbol
    unique.symbol <- subset(unique.symbol, end-start >= 300)
    dim(subset(unique.symbol, refseq == ""))

    if (!is.null(linc.file)) {
      linc <- readBed(linc.file)
      ## look at overlapping lincs and MGI symbol genes, and consolidate
      inter <- as.data.frame(findOverlaps(GRangesFromBed(linc), GRangesFromBed(unique.symbol[,1:6])))
      #unique.symbol[inter$subjectHits,]
      lincs.to.remove <- c()
      genes.to.remove <- c()
      to.remove <- mapply(function(lincHit, geneHit) {
        if (greplany(c("Gm","Rik"), unique.symbol$symbol[geneHit]))
          genes.to.remove <<- c(genes.to.remove, geneHit)
        else
          lincs.to.remove <<- c(lincs.to.remove, lincHit)
        linc$end[lincHit] - unique.symbol$start[geneHit]
        }, inter$queryHits, inter$subjectHits)

      linc.to.use <- linc[-lincs.to.remove,1:6]
      linc.to.use$symbol <- linc.to.use$refseq <- NA
      unique.symbol <- rbind(linc.to.use, unique.symbol[-genes.to.remove,])
    }

    genes <- unique.symbol[order(unique.symbol$start),1:8]
    genes <- genes[order(genes$chr),]
    if (write.file) write.table(genes, file=file, sep='\t', quote=F, col.names=T, row.names=F)
    return(genes)
  }
}


getESCAnnotatedGeneSet <- function(genome="mm9", file="/seq/lincRNA/data/mm9/mm9.ESCAnnotatedGeneData.txt", write.file=FALSE) {
  ## December 10, 2015
  ## including both mm9 protein-coding and lincRNA genes, one isoform per gene, just genomic span
  if (file.exists(file)) {
    genes <- read.delim(file, stringsAsFactors=F)
    return(genes)
  } else {
    gene.map <- getSimpleGeneNameMap()

    ## Annotate with expression level
    expression <- read.delim("/seq/lincRNA/RAP/PluripotencyExperiment/151002_FractionationFlavopiridol/analysis/exonRpkmTable.txt")
    gene.map$pA.WCE.RPKM <- mapply(function(name, refseq) {
      result <- numeric(0)
      if (refseq != "" & !is.na(refseq))
        result <- expression[grep(refseq, rownames(expression)),1]

      if (length(result) < 1) {
        result <- expression[grep(name, rownames(expression)), 1]
      }

      if (length(result) > 1) {
        cat(paste("Found more than one entry for ", name, refseq,"\n"))
        result <- result[1]
      }

      if (length(result) < 1) result <- NA
      return(result)
    }, gene.map$name, gene.map$refseq)


    ## Annotate with transcription level (non-polyA chromatin RPKM)
    transcription <- readCountReads("/seq/lincRNA/RAP/PluripotencyExperiment/151002_FractionationFlavopiridol/aligned/Total-Chromatin-DMSO/Total-Chromatin-DMSO.unique_alignment.CountReads.bed")
    gene.map$transcription.RPKM <- mapply(function(name, refseq) {
      result <- numeric(0)
      if (refseq != "" & !is.na(refseq))
        result <- transcription$rpkm[grep(refseq, transcription$name)]

      if (length(result) < 1) {
        result <- transcription$rpkm[grep(name, transcription$name)]
      }

      if (length(result) > 1) {
        cat(paste("Found more than one entry for ", name, refseq,"\n"))
        result <- result[1]
      }

      if (length(result) < 1) result <- NA
      return(result)
    }, gene.map$name, gene.map$refseq)

    ## Calculate TSS
    gene.map <- addTSSToBED(gene.map)

    ## Bivalent promoter in ESCs
    library(GenomicRanges)
    hmm <- readBed("/seq/lincRNA/RAP/Promoters/150827_TRIP_Analysis/mESC_cStates_HMM.sorted.bed")
    hmm.poised <- GRangesFromBed(subset(hmm, name == "10_Poised_Promoter")[,1:3])
    overlap <- findOverlaps(GRangesFromBed(with(gene.map, data.frame(chr, tss-1000, tss+1000))), hmm.poised)
    gene.map$bivalent <- FALSE
    gene.map$bivalent[queryHits(overlap)] <- TRUE

    ## Add H3K27ac total signal around promoter (2 kb on either side)
    to.write <- with(gene.map, data.frame(chr=chr, start=tss-2000, end=tss+2000, name=name))
    write.table(to.write, file="/seq/lincRNA/Jesse/tmp/tmp.H3K27ac.bed", sep='\t', quote=F, row.names=F, col.names=F)
    system(paste0("source /broad/software/scripts/useuse; reuse Java-1.7; ",
      "java -Xmx8g -cp /seq/lincRNA/Jesse/bin/scripts/Nextgen_140308.jar broad.pda.seq.rap.CountReads TARGET=/seq/lincRNA/RAP/PluripotencyExperiment/data/ChIP/SRR1202459.sorted.bam OUTPUT=/seq/lincRNA/RAP/PluripotencyExperiment/data/ChIP/SRR1202459.TSS.bed ANNOTATION_FILE=/seq/lincRNA/Jesse/tmp/tmp.H3K27ac.bed SCORE=count MASK_FILE=null PAIRED_END=false SIZES=/seq/lincRNA/data/mm9/sizes"
    ))
    result <- readCountReads("/seq/lincRNA/RAP/PluripotencyExperiment/data/ChIP/SRR1202459.TSS.bed")
    colnames(result)[8] <- "symbol"
    gene.map$H3K27ac.tss <- sapply(gene.map$name, function(x) {
      curr.name <- subset(result, name == x)
      if (nrow(curr.name) == 0) return(NA)
      else return(curr.name$score[1])
    })
    sum(is.na(gene.map$H3K27ac.tss))

    gene.map$expressed <- with(gene.map, (transcription.RPKM > 0.5) | (pA.WCE.RPKM > 0.5))

    final <- gene.map
    if (write.file) write.table(final, file=file, sep='\t', quote=F, col.names=T, row.names=F)
    return(final)
  }
}



getAnnotatedGeneSet <- function(cell.type="K562", write.file=FALSE) {
  ## April 22, 2016
  ## generalized version of above
  ## including both mm9 protein-coding and lincRNA genes, one isoform per gene, just genomic span

  if (cell.type == "mESC") {
    genome="mm9"
    file="/seq/lincRNA/data/mm9/mm9.ESCGeneSetWithLncRNA.txt"
    rpkm.file="/seq/lincRNA/RAP/PluripotencyExperiment/151002_FractionationFlavopiridol/analysis/exonRpkmTable.txt"
    txn.file="/seq/lincRNA/RAP/PluripotencyExperiment/151002_FractionationFlavopiridol/aligned/Total-Chromatin-DMSO/Total-Chromatin-DMSO.unique_alignment.CountReads.bed"
    hmm.file="/seq/lincRNA/RAP/Promoters/150827_TRIP_Analysis/mESC_cStates_HMM.sorted.bed"
    bivalent.state="10_Poised_Promoter"
    sizes="/seq/lincRNA/data/mm9/sizes"
    h3k27ac.bam="/seq/lincRNA/RAP/PluripotencyExperiment/data/ChIP/SRR1202459.sorted.unique_alignment.bam"
    txn.rpkm.cutoff=0.5
    WCE.pA.cutoff=0.5
  } else if (cell.type == "K562") {
    genome="hg19"
    file="/seq/lincRNA/data/hg19/hg19.K562AnnotatedGeneSet.txt"
    rpkm.file=NULL
    txn.file=NULL
    hmm.file="/seq/lincRNA/Jesse/CRISPR_Screen/data/wgEncodeBroadHmmK562HMM.bed"
    bivalent.state="3_Poised_Promoter"
    sizes="/seq/lincRNA/data/hg19/sizes"
    h3k27ac.bam="/seq/lincRNA/Jesse/CRISPR_Screen/data/wgEncodeBroadHistoneK562H3K27ac_ENCFF000BWZ.bam"
    txn.rpkm.cutoff=0.5
    WCE.pA.cutoff=0.5
  } else if (cell.type == "GM12878") {
    genome="hg19"
    file="/seq/lincRNA/data/hg19/hg19.GM12878AnnotatedGeneSet.txt"
    rpkm.file=NULL
    txn.file=NULL
    hmm.file="/seq/lincRNA/data/hg19/ENCODE/GM12878/wgEncodeBroadHmmGm12878HMM.bed"
    bivalent.state="3_Poised_Promoter"
    sizes="/seq/lincRNA/data/hg19/sizes"
    h3k27ac.bam="/seq/lincRNA/data/hg19/ENCODE/GM12878/wgEncodeBroadHistoneGm12878H3k27acStdAlnRep1.bam"
    txn.rpkm.cutoff=0.5
    WCE.pA.cutoff=0.5
  } else if (cell.type == "Karpas422") {
    genome="hg19"
    file="/seq/lincRNA/data/hg19/hg19.Karpas422AnnotatedGeneSet.txt"
    rpkm.file=NULL
    txn.file=NULL
    hmm.file=NULL
    bivalent.state="3_Poised_Promoter"
    sizes="/seq/lincRNA/data/hg19/sizes"
    h3k27ac.bam="/seq/lincRNA/data/hg19/ENCODE/Karpas422/Karpas422-H3K27ac-Alignment_Post_Processing_7601.bam"
    txn.rpkm.cutoff=0.5
    WCE.pA.cutoff=0.5
  } else if (cell.type == "HUVEC") {
      genome="hg19"
    file="/seq/lincRNA/data/hg19/hg19.HUVECAnnotatedGeneSet.txt"
    rpkm.file=NULL
    txn.file=NULL
    hmm.file="/seq/lincRNA/data/hg19/ENCODE/HUVEC/wgEncodeBroadHmmHuvecHMM.bed"
    bivalent.state="3_Poised_Promoter"
    sizes="/seq/lincRNA/data/hg19/sizes"
    h3k27ac.bam="/seq/lincRNA/data/hg19/ENCODE/HUVEC/wgEncodeBroadHistoneHuvecH3k27acStdAlnRep1.bam"
    txn.rpkm.cutoff=0.5
    WCE.pA.cutoff=0.5
  } else {
    stop(paste0("Genome ",genome," not recognized"))
  }

  if (file.exists(file)) {
    genes <- read.delim(file, stringsAsFactors=F)
    return(genes)
  } else {
    gene.map <- getSimpleGeneNameMap(genome=genome)

    ## Annotate with expression level
    if (!is.null(rpkm.file)) {
      expression <- read.delim(rpkm.file)
      gene.map$pA.WCE.RPKM <- mapply(function(name, refseq) {
        result <- numeric(0)
        if (refseq != "" & !is.na(refseq))
          result <- expression[grep(refseq, rownames(expression)),1]

        if (length(result) < 1) {
          result <- expression[grep(name, rownames(expression)), 1]
        }

        if (length(result) > 1) {
          cat(paste("Found more than one entry for ", name, refseq,"\n"))
          result <- result[1]
        }

        if (length(result) < 1) result <- NA
        return(result)
      }, gene.map$name, gene.map$refseq)
    } else gene.map$pA.WCE.RPKM <- NA

    ## Annotate with transcription level (non-polyA chromatin RPKM)
    if (!is.null(txn.file)) {
      transcription <- readCountReads(txn.file)
      gene.map$transcription.RPKM <- mapply(function(name, refseq) {
        result <- numeric(0)
        if (refseq != "" & !is.na(refseq))
          result <- transcription$rpkm[grep(refseq, transcription$name)]

        if (length(result) < 1) {
          result <- transcription$rpkm[grep(name, transcription$name)]
        }

        if (length(result) > 1) {
          cat(paste("Found more than one entry for ", name, refseq,"\n"))
          result <- result[1]
        }

        if (length(result) < 1) result <- NA
        return(result)
      }, gene.map$name, gene.map$refseq)
    } else gene.map$transcription.RPKM <- NA

    ## Calculate TSS
    gene.map <- addTSSToBED(gene.map)

    ## Bivalent promoter
    if (!is.null(hmm.file)) {
      library(GenomicRanges)
      hmm <- readBed(hmm.file)
      hmm.poised <- GRangesFromBed(subset(hmm, name == bivalent.state)[,1:3])
      overlap <- findOverlaps(GRangesFromBed(with(gene.map, data.frame(chr, tss-1000, tss+1000))), hmm.poised)
      gene.map$bivalent <- FALSE
      gene.map$bivalent[queryHits(overlap)] <- TRUE
    }

    ## Add H3K27ac total signal around promoter (2 kb on either side)
    to.write <- with(gene.map, data.frame(chr=chr, start=tss-2000, end=tss+2000, name=name))
    write.table(to.write, file="/seq/lincRNA/Jesse/tmp/tmp.H3K27ac.bed", sep='\t', quote=F, row.names=F, col.names=F)
    system(paste0("source /broad/software/scripts/useuse; reuse Java-1.7; ",
      "java -Xmx8g -cp /seq/lincRNA/Jesse/bin/scripts/Nextgen_140308.jar broad.pda.seq.rap.CountReads TARGET=",
      h3k27ac.bam," OUTPUT=/seq/lincRNA/Jesse/tmp/tmp.H3K27ac.CountReads.bed VALIDATION_STRINGENCY=LENIENT ANNOTATION_FILE=/seq/lincRNA/Jesse/tmp/tmp.H3K27ac.bed SCORE=count MASK_FILE=null PAIRED_END=false SIZES=/seq/lincRNA/data/mm9/sizes"
    ))
    result <- readCountReads("/seq/lincRNA/Jesse/tmp/tmp.H3K27ac.CountReads.bed")
    gene.map$H3K27ac.tss <- sapply(gene.map$name, function(x) {
      curr.name <- subset(result, name == x)
      if (nrow(curr.name) == 0) return(NA)
      else return(curr.name$score[1])
    })
    sum(is.na(gene.map$H3K27ac.tss))

    if (!is.null(txn.file) & !is.null(rpkm.file)) {
      gene.map$expressed <- with(gene.map, (transcription.RPKM > txn.rpkm.cutoff) | (pA.WCE.RPKM > pA.WCE.cutoff))
    } else {
      gene.map$expressed <- NA
    }

    final <- gene.map
    if (write.file) write.table(final, file=file, sep='\t', quote=F, col.names=T, row.names=F)
    return(final)
  }
}


addTSSToBED <- function(bed) {
  bed$tss <- bed$start
  bed$tss[bed$strand == "-"] <- bed$end[bed$strand == "-"]
  ## Ideally would grab the TSS from CRISPRi defined files
  bed$tss[bed$name == "Cd276"] <- 58402943 ## this is incorrectly called
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

###########################################################################################
## Make a count table from CountReads files

makeCountTable <- function(sample.names, aligned.dir, tail=".unique_alignment.CountReads.bed", col="score", no.sample.dir=FALSE) {
  sample.names <- as.character(as.matrix(sample.names))
  countTable <- c()
  cols <- c()
  for (i in 1:length(sample.names)) {
    tryCatch({
      sample <- sample.names[i]
      if (no.sample.dir) 
        curr.file <- paste0(aligned.dir,"/",sample,tail)
      else
        curr.file <- paste0(aligned.dir,"/",sample,"/",sample,tail)
      curr <- readCountReads(curr.file)
      countTable <- cbind(countTable, with(curr, get(col)))
      rownames(countTable) <- as.character(as.matrix(curr$name))
      cols <- c(cols, sample)
    }, error = function(e) { print(e); traceback() }
          )
  }
  colnames(countTable) <- cols
  return(countTable)
}


ternaryPlotJuicer <- function(x, ...) {
  ## Jesse's version of ternaryplot -- plot just calls ternary plot, but this function also
  ## returns the coordinates of the points for additional plotting purposes
  require(vcd)
  ternaryplot(x, ...)
  s <- rowSums(x)
  x <- x/s
  top <- sqrt(3)/2
  xp <- x[, 2] + x[, 3]/2
  yp <- x[, 3] * top
  return(cbind(xp, yp))
}


barplotWithErrorBars <- function(data, bar.col, error.col, name.col, ...) {
  h <- barplot(data[,bar.col], ...)
  axis(1, labels=data[,name.col], at=h, las=2)
  abline(h=1, lty=2, col='gray')
  #abline(v=0, lty=1, col='black')
  segments(h, data[,bar.col] - data[,error.col], h, data[,bar.col] + data[,error.col])
}


barplotGroupsWithErrorBars <- function(mat, error.mat=NULL, ...) {
    h <- barplot(mat, beside=TRUE, ...)
    if (!is.null(error.mat)) segments(h, mat - error.mat, h, mat + error.mat)
    return(h)
}


geometricMean  <- function(x, na.rm=TRUE) {
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


fixPythonLogicalFields <- function(df) {
  for (i in 1:ncol(df)) {
    if (all(df[,i] %in% c("True","False",""))) {
      df[,i] <- as.logical(as.character(as.matrix(df[,i])))
    }
  }
  return(df)
}