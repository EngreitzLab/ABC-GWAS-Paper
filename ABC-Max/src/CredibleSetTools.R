## Jesse Engreitz
## December 9, 2018
## Functions to support GWAS - EP prediction analysis

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))


## for dplyr < 0.7
## TODO: Fix conda env to use later version of R and dplyr
pull <- function(x,y) {x[,if(is.name(substitute(y))) deparse(substitute(y)) else y, drop = FALSE][[1]]}

####################################################################################
## Functions to load in variant predictions from ABC

loadVariantOverlap <- function(
  overlap.file, 
  genes.uniq, 
  genes, 
  variant.names=NULL, 
  overwriteTSS=FALSE, 
  isTargetGeneTSS=TRUE) {

  ## Loads overlap file, and filters to variants in variant.names
  x <- read.delim(gzfile(overlap.file), check.names=F)
  
  if (!is.null(variant.names)) {
    tmp <- data.frame(variant=factor(as.character(as.matrix(variant.names)), levels=levels(x$QueryRegionName)))
    tmp <- subset(tmp, !is.na(variant))
    x <- merge(x, tmp, by.x="QueryRegionName", by.y="variant")
  }
  
  ## Merge and add various other annotations
  if (overwriteTSS) {
    if ("TargetGeneTSS" %in% colnames(x)) x <- x %>% select(-TargetGeneTSS)
    x <- merge(x, with(genes.uniq, data.frame(TargetGene=name, TargetGeneTSS=tss)), by="TargetGene")
  }
  
  if (isTargetGeneTSS) {
    x$isOwnTSS <- with(x, TargetGeneTSS >= start & TargetGeneTSS <= end)
  }
  codingSymbols <- subset(genes, grepl(";NM_", name))$symbol
  
  x$TargetGeneIsCoding <- as.character(as.matrix(x$TargetGene)) %in% codingSymbols
  return(x)
}


getCodingGenes <- function(genes=NULL, genes.uniq=NULL) {
  if (is.null(genes)) {
    genes <- readBed("/seq/lincRNA/data/hg19/RefSeqCurated.170308.bed")
    genes$symbol <- unlist(lapply(strsplit(as.character(as.matrix(genes$name)), ";"), "[", 1))
  }
  if (is.null(genes.uniq)) genes.uniq <- readBed("/seq/lincRNA/data/hg19/RefSeqCurated.170308.bed.CollapsedGeneBounds.bed")
  codingSymbols <- as.character(as.matrix(subset(genes, grepl(";NM_", name))$symbol))
  result <- subset(genes.uniq, name %in% codingSymbols)
  result$tss <- result$start; result$tss[result$strand == "-"] <- result$end[result$strand == "-"]
  return(result)
}


filterVariantOverlap <- function(overlap, predCol, cutoff, tss.cutoff, ignore.list=c()) {
  ## Implements a different cutoff for distal enhancers versus distal promoters
  if (!is.na(cutoff)) 
    overlap <- subset(overlap, ((class != "tss" & class != "promoter") | isOwnTSS) | (get(predCol) >= tss.cutoff))
  if (!is.na(tss.cutoff))
    overlap <- subset(overlap, ((class == "tss" | class == "promoter")) | (get(predCol) >= cutoff))

  overlap <- overlap %>% filter( !(TargetGene %in% ignore.list) )
  
  return(overlap)
}


annotateVariantOverlaps <- function(overlap, variant.list, all.cs, var.cols=c("CredibleSet","Disease","PosteriorProb","Coding","SpliceSite","Promoter","LocusID")) {
  cols.to.remove <- which(colnames(variant.list) %in% c("chr","position","start","end"))  
  all.flat <- overlap
  all.flat <- merge(all.flat, variant.list[,-cols.to.remove], by.x="QueryRegionName", by.y="variant")
  return(all.flat)
}



#######################################################################################
## Functions for making and loading data from permutation tables written by MakeVariantCountTables.R

getVariantByCellsTable <- function(overlap, score.col="PosteriorProb", isTargetGene=TRUE, isCellType=TRUE) {
    if ( isCellType & isTargetGene ) {
  	variant.by.cells <- overlap %>% group_by(QueryRegionName,CellType) %>% summarise( n.genes=n(), max.ABC=max(ABC.Score), TargetGenes=paste(TargetGene,collapse=','), PosteriorProb=max(PosteriorProb) ) %>% as.data.frame()
	return(variant.by.cells)
  } else {
	  if ( !isTargetGene & !isCellType ) {
  		variant.by.cells <- overlap %>% group_by(QueryRegionName) %>% summarise( n.genes=n(), PosteriorProb=max(PosteriorProb) ) %>% as.data.frame()
	  } else if (isCellType) {
		  variant.by.cells <- overlap %>% group_by(QueryRegionName,CellType) %>% summarise( n.genes=n(), PosteriorProb=max(PosteriorProb) ) %>% as.data.frame()
	  } else if (isTargetGene) {
		  variant.by.cells <- overlap %>% group_by(QueryRegionName) %>% summarise( n.genes=n(), TargetGenes=paste(TargetGene,collapse=','), PosteriorProb=max(PosteriorProb) ) %>% as.data.frame()}
	}
	return(variant.by.cells)
}


countCellTypeOverlaps <- function(dat, cell.type.list, variant.names=NULL, cell.groups=NULL, weightByPIP=FALSE) {
  cell.type.list <- na.omit(cell.type.list)
  
  if (!is.null(variant.names)) {
    dat <- dat %>% filter(QueryRegionName %in% variant.names)
  }
  
  if (nrow(dat) == 0) {
    result <- data.frame(CellType=cell.type.list, n=0)
  } else {
    
    if (!is.null(cell.groups)) {
      ## Convert cell types to cell type groups
      df <- data.frame(CellType=cell.type.list, CellGroup=cell.groups)
      dat <- merge(dat, df, by="CellType")
      dat$CellType <- dat$CellGroup
      cell.type.list <- unique(cell.groups)
    }
    
    dat$CellType <- refactor(dat$CellType, cell.type.list)
    
    if (!weightByPIP) {
      dat <- unique(dat[,c("QueryRegionName","CellType")])
      result <- dat %>% group_by(CellType) %>% tally()
    } else {
      dat <- unique(dat[,c("QueryRegionName","CellType","PosteriorProb")])
      result <- dat %>% group_by(CellType) %>% summarise(n=sum(PosteriorProb))
    }
    result <- merge(data.frame(CellType=cell.type.list), result, all.x=T)
    result$n[is.na(result$n)] <- 0
  }
  return(result)
}


computeCellTypeEnrichment <- function(variants.by.cells, variant.list, cell.type.annot, trait, cell.group.by=NULL, score.col="PosteriorProb", min.score=0.1, bg.vars=NULL, bg.overlap=NULL, noPromoter=FALSE, isCellType=TRUE, enrichment.threshold=0.001) {
  ## Computes the enrichment of variants with high vs low posterior probabilities in each cell type
  ## variants.by.cells    data.frame output by getVariantByCellsTable
  ## variant.list         data.frame with all variant info. Will subset the variants.by.cells df by this list of variants
  ## cell.type.annot      data.frame containing cell types and columns with categorical annotations
  ## cell.group.by        column name from cell.type.annot to group the enrichment calculation
  
  # If a score threshold is probivided, selecting significant variants. Else, using all variants.
  if (!(is.null(min.score))){
    hi.vars <- subset(variant.list, get(score.col) >= min.score)$variant
  } else {
    hi.vars <- variant.list$variant
  }
  
  ##TODO: Rewrite this using interpretable variable names
  bg_V4 <- bg.vars$V4
  
  stopifnot(nrow(hi.vars) > 0)
  stopifnot(nrow(bg.vars) > 0)
  
  message(cell.group.by)
  group.list <- if (!is.null(cell.group.by)) cell.type.annot[[cell.group.by]] else NULL
  if (isCellType){
  	hi.count <- countCellTypeOverlaps(variants.by.cells, cell.type.list=cell.type.annot, variant.names=hi.vars, cell.groups=group.list)
  	bg_V8_count <- bg.overlap %>% distinct(V4,V8)
	bg_V8 <- bg_V8_count %>% group_by(V8) %>% tally()
	hi.count$n.ctrl <- bg_V8$n
  } else {
	  hi.count  <- variants.by.cells %>% filter(QueryRegionName %in% hi.vars)
	  hi.count <- data.frame(n=length(hi.count$QueryRegionName))
	  bg_V8 <- bg.overlap %>% distinct(V4)
	  hi.count$n.ctrl <- length(unique(bg_V8$V4))
  }

  hi.count$total <- length(hi.vars)
  hi.count$total.ctrl <- length(bg.vars$V4)
  hi.count$Disease <- trait
  hi.count$prop.snps <- with(hi.count, n/total)
  hi.count$enrichment <- with(hi.count, n/total / ( (n.ctrl)/total.ctrl ))
  hi.count$log10.p <- -with(hi.count, mapply(FUN=phyper, n, total, total.ctrl, n+n.ctrl, log.p=TRUE, lower.tail=F))/log(10)
  hi.count$log10.p[is.infinite(hi.count$log10.p)] <- NA
  hi.count$p <- with(hi.count, mapply(FUN=phyper, n, total, total.ctrl, n+n.ctrl, log.p=FALSE, lower.tail=F))
  hi.count$Significant <- p.adjust((hi.count$p), method="bonferroni") <= enrichment.threshold
  hi.count$enrichment[with(hi.count, n==0 & n.ctrl==0)] <- NA  
  return(hi.count)
}


####################################################################################
## Gene Prioritization Table

getGenePrioritizationTable <- function(
  all.clean, 
  all.cs, 
  genes, 
  genes.uniq, 
  cell.types, 
  cell.type.annot, 
  score.col="ABC.Score",
  score.min=0.015,
  var.score.col="PosteriorProb",
  var.score.min=0.1,
  max.distance=1000000,
  method.name="ABC") { 


  cs.tables <- list()
  for (cs.name in all.cs$CredibleSet) {
    cs.tables[[cs.name]] <- getGenePrioritizationTableForCredibleSet(
        cs.name, 
        all.clean, 
        all.cs,
        genes, 
        genes.uniq, 
        cell.types, 
        cell.type.annot, 
        score.col,
        score.min,
        var.score.col,
        var.score.min,
        max.distance,
        method.name)
  }
  result <- do.call(rbind, cs.tables) 
  return(result)
}


writeGenePrioritizationTable <- function(gp, file) {
  write.tab(
'# Column definitions:   
# trait : GWAS phenotype
# CredibleSet : credible set ID for trait
# CredibleSet.AnyCodingVariant : TRUE if credible set contains a variant overlapping a coding sequence
# CredibleSet.AnySpliceSiteVariant : TRUE if credible set contains a variant within 10bp of a splice site of a coding gene
# CredibleSet.NoncodingWithSigVariant : TRUE if credible set does not contain a coding or splice site variant, and if it contains at least one variant with variant score above the specified threshold
# CredibleSet.AnyPrediction : TRUE if credible set has a prediction for this method
# CredibleSet.nNearbyGenes : Number of genes within distance threshold
# TargetGene : Gene Symbol
# TargetGene.chr : chromosome
# TargetGene.TSS : coordinate of gene transcription start site
# DistanceToTSS : Distance from best variant to the transcription start site
# DistanceToTSSRank : Rank among genes of the distance from best variant to the transcription start site
# Distance : Distance from best variant to the gene body
# DistanceRank : Rank among genes of the distance from best variant to the gene body
# GeneScore.* : Max score for this gene across cell types and variants
# CellTypes.* : Cell types in which the prediction is made
# Variants.* : Variants with variant scores above threshold responsible for the prediction
# GeneRank.* : Rank of the prediction score for this gene among other genes in the credible set
# GenePrediction.* : TRUE if gene is the gene with the max by this predictor for the credible set
# GenePredictionMax.* : TRUE if gene has the maximum score for this credible set',

  file=file, col.names=F)
  write.tab(gp, file=file, append=TRUE)
}

getGenePrioritizationTableForCredibleSet <- function(
  cs.name, 
  all.clean, 
  all.cs,
  genes, 
  genes.uniq, 
  cell.types, 
  cell.type.annot, 
  score.col,
  score.min,
  var.score.col,
  var.score.min,
  max.distance, 
  method.name) { 

  curr.cs <- subset(all.cs, CredibleSet == cs.name)
  stopifnot(nrow(curr.cs) == 1)

  cs.genes <- getGenesNearCredibleSet(
    genes.uniq, 
    all.clean, 
    score.col, 
    score.min, 
    var.score.col, 
    var.score.min, 
    cs.name, 
    curr.cs$chr, 
    curr.cs$start, 
    curr.cs$end, 
    max.distance)

  stopifnot(length(unique(all.cs$Disease)) == 1)

  gp <- subset(genes.uniq, name %in% cs.genes)
  gp <- gp %>% transmute(
    trait=unique(all.cs$Disease),
    CredibleSet=cs.name,
    CredibleSet.AnyCodingVariant=curr.cs$AnyCoding,
    CredibleSet.AnySpliceSiteVariant=curr.cs$AnySpliceSite,
    CredibleSet.NoncodingWithSigVariant=with(curr.cs, !AnyCoding & !AnySpliceSite & MaxVariantScore >= var.score.min),
    CredibleSet.AnyPrediction=NA,
    CredibleSet.nNearbyGenes=n(),
    TargetGene=name,
    TargetGene.chr=chr,
    TargetGene.TSS=tss,
    DistanceToTSS=abs(curr.cs$BestSNPPos-tss),
    DistanceToTSSRank=rank(DistanceToTSS, ties.method="min"), 
    Distance=distance(IRanges(curr.cs$BestSNPPos,curr.cs$BestSNPPos),IRanges(start,end-1)),
    DistanceRank=rank(Distance, ties.method="min"))

  curr.pred <- subset(all.clean, CredibleSet == cs.name)
  curr.pred <- subset(curr.pred, CellType %in% cell.types)
  curr.pred <- subset(curr.pred, TargetGene %in% genes.uniq$name)
  if (!is.null(score.col) & nrow(curr.pred) > 0) curr.pred <- subset(curr.pred, get(score.col) >= score.min)
  if (!is.null(var.score.col) & nrow(curr.pred) > 0) curr.pred <- subset(curr.pred, get(var.score.col) >= var.score.min)

  if (nrow(curr.pred) > 0) {
    pred.scores <- do.call(rbind, 
      by(curr.pred, curr.pred$TargetGene, function(x) {
        df <- data.frame(TargetGene=x$TargetGene[1])
        df[[paste0("GeneScore.",method.name)]] <- max(x[,score.col])
        df[[paste0("CellTypes.",method.name)]] <- paste0(sort(unique(x$CellType)),collapse=',')
        df[[paste0("Variants.",method.name)]] <- paste0(sort(unique(x$QueryRegionName)),collapse=',')
        return(df)
      }, simplify=FALSE))

    gp.scores <- merge(gp, pred.scores, by="TargetGene", all.x=TRUE)[,union(names(gp), names(pred.scores))]
    gp.scores[[paste0("GeneRank.",method.name)]] <- rank(-gp.scores$GeneScore, ties.method="min")
    gp.scores[[paste0("GenePrediction.",method.name)]] <- !is.na(gp.scores[,paste0("GeneScore.",method.name)])
    gp.scores[[paste0("GenePredictionMax.",method.name)]] <- !is.na(gp.scores[,paste0("GeneScore.",method.name)]) & (gp.scores[,paste0("GeneRank.",method.name)] == 1 | is.na(gp.scores[,paste0("GeneScore.",method.name)]))

  } else {
    gp.scores <- gp
    gp.scores[[paste0("GeneScore.",method.name)]] <- NA
    gp.scores[[paste0("CellTypes.",method.name)]] <- NA
    gp.scores[[paste0("Variants.",method.name)]] <- NA
    gp.scores[[paste0("GeneRank.",method.name)]] <- NA
    gp.scores[[paste0("GenePrediction.",method.name)]] <- FALSE
    gp.scores[[paste0("GenePredictionMax.",method.name)]] <- FALSE
  }
  gp.scores <- gp.scores %>% arrange(TargetGene.TSS)
  return(gp.scores)
}

getGenesNearCredibleSet <- function(
  genes.uniq, 
  all.flat, 
  score.col, 
  score.min, 
  var.score.col, 
  var.score.min, 
  cs.name, 
  cs.chr, 
  cs.start, 
  cs.end, 
  max.distance) {

  nearby.genes <- genes.uniq %>% filter(chr == cs.chr & tss >= cs.start-max.distance & tss <= cs.end+max.distance) %>% pull(name) %>% as.matrix() %>% as.character()
  abc.genes <- all.flat %>% filter(CredibleSet == cs.name)
  if (!is.null(score.col) & nrow(abc.genes) > 0) 
    abc.genes <- subset(abc.genes, get(score.col) >= score.min)
  if (!is.null(var.score.col) & nrow(abc.genes) > 0) 
    abc.genes <- subset(abc.genes, get(var.score.col) >= var.score.min)
  abc.genes <- abc.genes$TargetGene %>% as.matrix() %>% as.character()

  cs.genes <- unique(c(nearby.genes, abc.genes))
  return(cs.genes)
}


getPredictedCellTypes <- function(gp, cellTypes) {
  apply(gp, 1, function(row) {
    scores <- as.numeric(row[cellTypes]); names(scores) <- cellTypes
    scores <- scores[scores > 0]
    scores <- sort(scores)
    paste0(names(scores), collapse=',')
  })
}

getPredictedVariants <- function(gp, all.flat, cellTypes, score.col="PosteriorProb", min.score=0.1) {
  tmp <- all.flat %>% filter(CellType %in% cellTypes & get(score.col) >= min.score) %>% group_by(CredibleSet, TargetGene) %>% summarise(Variants=paste0(unique(QueryRegionName), collapse=',')) %>% as.data.frame()
  tmp <- unfactor(tmp)
  gp$origOrder <- 1:nrow(gp)
  gp.tmp <- merge(gp, tmp, all.x=TRUE, by=c("CredibleSet","TargetGene"))
  return(gp.tmp$Variants[order(gp.tmp$origOrder)])
}

# cellTypeFlag="AnyDisease_FMOverlap_Enriched" aka "IBD_FMOverlap_Enriched"
# TODO: get correct cell type column from main script
getABCMaxTable <- function(gp, all.flat, cell.type.annot, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold, cellTypeFlag="IBD_FMOverlap_Enriched") {
  gp.pred <- subset(gp, get(paste0("ConnectionStrengthRank.Binary.",cellTypeFlag)) == 1)
  enriched.cell.types <- as.character(as.matrix(subset(cell.type.annot, get(paste0("Binary.",cellTypeFlag)))$CellType))
  gp.pred$CellTypes <- getPredictedCellTypes(gp.pred, enriched.cell.types)
  gp.pred$Variants <- getPredictedVariants(gp.pred, all.flat, enriched.cell.types, score.col=opt$variantScoreCol, min.score=opt$variantScoreThreshold)
  abcmax <- gp.pred[,c("Disease","TargetGene","CredibleSet","Variants","CellTypes")]
  abcmax <- abcmax[order(abcmax$TargetGene),]
  return(abcmax)
}


checkInteger <- function(n) n == round(n)


getPrecisionRecall <- function(pred, labels, weights=NULL, ...) {
  ## Used currently by CollateABCTrainingData.R
  stopifnot(is.logical(labels))
  stopifnot(length(pred) == length(labels))
  
  if (is.null(weights)) weights <- rep(1, length(pred))
  
  if (is.logical(pred)) {
    recall <- sum(weights[pred & labels], na.rm=T)/sum(weights[labels], na.rm=T)
    precision <- sum(weights[pred & labels], na.rm=T)/sum(weights[pred], na.rm=T)
    result <- data.frame(
      Precision=precision,
      Recall=recall,
      nPred=sum(weights[pred],na.rm=T),
      nTruePred=sum(weights[pred & labels],na.rm=T),
      nTrue=sum(weights[labels],na.rm=T),
      alpha=TRUE,
      ...
    )
  } else if (all(checkInteger(pred), na.rm=T)) {
    ## Assume this is a ranking and that lower is better
    result <- do.call(rbind, lapply(na.omit(sort(unique(pred))), function(cutoff) {
      tmp <- getPrecisionRecall(pred <= cutoff, labels, weights, ...)
      tmp$alpha <- cutoff
      return(tmp)
    }))
  } else if (is.numeric(pred)) {
    result <- do.call(rbind, lapply(sort(unique(pred)), function(cutoff) {
      tmp <- getPrecisionRecall(pred >= cutoff, labels, weights, ...)
      tmp$alpha <- cutoff
      return(tmp)
    }))
  } else {
    stop("Class of 'pred' not supported: ", class(pred))
  }
  return(result)
}





