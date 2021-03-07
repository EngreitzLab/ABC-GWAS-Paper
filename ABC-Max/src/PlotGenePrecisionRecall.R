# Plot precision-recall at identifying known IBD genes


suppressPackageStartupMessages(library(optparse))


option_list <- list(
        make_option("--genePredTable", help="Input gene predictions file (for one or more traits)"),
        make_option("--outPdf", type="character", help="Output PDF for precision-recall plot"),
        make_option("--knownGenes", type="character", help="File with columns corresponding to lists of genes with which to evaluate predictors. Columns to use must start with 'GeneList.'"),
        make_option("--codeDir", type="character", default="ABC-Max-pipeline/", help="code directory"),
        make_option("--knownGeneMaxDistance", type="numeric", default=1000000, help="Maximum distance from credible set to TSS to count a known gene")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(RColorBrewer))

#######################################################################
## Check inputs

checkInputs <- function(opt) {
  if (!file.exists(opt$genePredTable)) stop(paste0("Gene prediction file does not exist: ", opt$genePredictionTable))
  if (!file.exists(opt$knownGenes)) stop(paste0("Known genes file does not exist: ", opt$knownGenes))
  if (!file.exists(opt$codeDir)) stop(paste0("Pipeline code directory does not exist: ", opt$codeDir))
}
checkInputs(opt)

source(paste0(opt$codeDir, "/src/Utilities.R"))
source(paste0(opt$codeDir, "/src/CredibleSetTools.R"))


#######################################################################
## Load gene prediction files 

gp <- read.delim(opt$genePredTable, check.names=F, stringsAsFactors=F, comment.char='#')
knownGenes <- read.delim(opt$knownGenes, check.names=F, stringsAsFactors=F, comment.char='#')

predictors <- colnames(gp)[greplany(c("GeneScore.","GenePrediction.","GenePredictionMax."), colnames(gp))]



#######################################################################
## Functions for plotting PR curves

mytheme <- theme_classic() + theme(
  axis.text = element_text(size = 13), 
  axis.title = element_text(size = 15))

getPrecisionBaseline <- function(gp) mean(1 / unique(gp[,c("CredibleSet","CredibleSet.nNearbyGenes")])$TotalNearbyGenes)

getPRTable <- function(gp, pred.cols) {
  pr <- do.call(rbind, c(
    with(gp, list(
      getPrecisionRecall(DistanceRank == 1, knownGene, Method="Closest Gene"),
      getPrecisionRecall(DistanceToTSSRank == 1, knownGene, Method="Closest TSS"))),
    lapply(pred.cols, function(col) getPrecisionRecall(gp[,col], gp$knownGene, Method=gsub("^Gene","",col)))
  ))
  return(pr)
}

getPRPlot <- function(pr, baseline=NULL, xlab="Recall") {
  pr <- pr %>% group_by(Method) %>% mutate(n=n())
  pr.points <- pr %>% filter(n==1)
  pr.lines <- pr %>% filter(n>1)
  p <- ggplot(pr.points, aes(x=Recall, y=Precision, color=Method)) 
  if (!is.null(baseline)) p <- p + geom_hline(yintercept=baseline, color='gray', linetype='dashed') 
  p <- p + geom_point(size=3)
  p <- p + geom_line(data=pr.lines)
  p <- p + mytheme + coord_fixed() + xlim(0,1) + ylim(0,1) + xlab(xlab) + ylab("Precision")
  return(p)
}

doOneKnownGeneList <- function(gene.list.name, gp, predictors, maxKnownGenes=1) {
  ## Current logic: Filter the gene prediction table to those credible sets with exactly one known gene nearby
  gp.plot <- gp %>% filter(DistanceToTSS <= opt$knownGeneMaxDistance) %>%
        mutate(knownGene=TargetGene %in% knownGenes[[gene.list.name]]) %>%
        group_by(CredibleSet) %>% mutate(nKnownGenes=sum(knownGene)) %>% ungroup() %>%
        filter(nKnownGenes > 0 & nKnownGenes <= maxKnownGenes & CredibleSet.NoncodingWithSigVariant) %>% as.data.frame()
  pr <- getPRTable(gp.plot, predictors)
  p <- getPRPlot(pr, baseline=getPrecisionBaseline(gp), xlab=paste0("Recall (n=",sum(gp.plot$knownGene),")"))
  print(p)
}


####################################################################
## PLOT:

pdf(file=opt$outPdf, width=6, height=4, onefile=T)
for (gl in colnames(knownGenes))
  doOneKnownGeneList(gl, gp, predictors)
dev.off()



