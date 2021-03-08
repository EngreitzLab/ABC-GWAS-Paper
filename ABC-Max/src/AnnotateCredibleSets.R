##################################################################
## Jesse Engreitz
## Annotate GWAS credible sets with enhancer-gene predictions

##  This executable script takes in a set of variant predictions (variants overlapped with ABC),
##  plus information about a set (or subset) of these variants and their credible sets,
##  and outputs a series of useful data tables that can be processed or analyzed in various
##  ways by other scripts.

suppressPackageStartupMessages(library(optparse))

# Parse options from command line
option.list <- list(
  make_option("--variants", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.variant.list.txt", help="File containing variants to consider"),
  make_option("--credibleSets", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/Huang2017-IBD/CredibleSets/IBDCombined.set1-2.cs.txt", help="File containing credible set annotations"),
  make_option("--predictionFile", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/abc.tsv.gz", help="Prediction overlap file (.txt.gz) output from intersecting variants with predicted enhancers"),
  make_option("--methodName", type="character", default="ABC", help="Name of E-G prediction method"),
  make_option("--predScoreCol", type="character", default="ABC.Score", help="Name of the column with the prediction score, e.g. ABC-score."),
  make_option("--outbase", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_out/IBD/", help="Output file basename"),
  make_option("--outEnrichment", type="character", default="./Enrichment.CellType.vsScore.tsv", help="Output cell type enrichment table"),
  make_option("--outGenePredTable", type="character", default="./GenePredictions.allCredibleSets.tsv", help="Output gene prediction table filename"),
  make_option("--cellType", type="character", default=TRUE, help="Do predictions have an associated cellType column?"),
  make_option("--TargetGene", type="character", default=TRUE, help="Do predictions have an associated targetGene column?"),
  make_option("--TargetGeneTSS", type="character", default=TRUE, help="Do predictions have an associated targetGeneTSS column?"),
  make_option("--minPredScore", type="numeric", default=NA, help="Cutoff on prediction score for distal elements"),
  make_option("--minPredScorePromoters", type="numeric", default=NA, help="Cutoff on prediction score for tss/promoter elements"),
  make_option("--trait", type="character", help="Name of the trait or disease"),
  make_option("--variantScoreCol", type="character", default="PosteriorProb", help="Score determining variant significance. If set to NULL, all variants in the input file will be included"),
  make_option("--variantScoreThreshold", type="numeric", default=0.1, help="Score cutoff for desired variants to analyze, e.g. PP>=0.1"),
  make_option("--backgroundVariants", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/all.bg.SNPs.bed.gz", help="A set of background variants to use in the enrichment analysis. Bed format with chr, start, end, rsID"),
  make_option("--bgOverlap", type="character", default="/oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_out/ABC.OverlapAllSNPs.tsv.gz", help="Background variant overlap with predictions"),
  make_option("--genes", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/RefSeqCurated.170308.bed", help="RefSeq gene BED file; this is to pull RefSeq IDs to determine coding/noncoding"),
  make_option("--genesUniq", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-MAX-pipeline/Test_data/RefSeqCurated.170308.bed.CollapsedGeneBounds.bed", help="Collapsed RefSeq gene BED file used for E-G predictions"),
  make_option("--genePredMaxDistance", type="numeric", default=1000000, help="Gene prediction table: Include genes within this distance"),
  make_option("--biosampleEnrichThreshold", type="numeric", default=0.001, help="Gene prediction and enrichment tables: Bonferroni-adjusted p-value to call biosample as significantly enriched for overlapping variants (set to 1 to include all biosamples in gene prediction table)"),

  make_option("--cellTypeTable", type="character", help="Table with annotations of cell types, with columns 'CellType', 'Categorical.*', 'Binary.*' for plotting enrichments"),
  make_option("--removePromoterVariants", type="logical", default=FALSE, help="Remove credible sets with promoter variants from the filter.cs list"),
  make_option("--removeCodingVariants", type="logical", default=TRUE, help="Remove credible sets with coding variants from the filter.cs list"),
  make_option("--removeSpliceSiteVariants", type="logical", default=TRUE, help="Remove credible sets with splice site variants from the filter.cs list"),
  make_option("--codeDir", type="character", default="//oak/stanford/groups/akundaje/kmualim/ABC-Max-pipeline/Utilities/", help="Directory to code base")
)

opt <- parse_args(OptionParser(option_list=option.list))
dput(opt)

# Loading libraries and utilities
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(source(paste0(opt$codeDir, "Utilities.R")))
suppressPackageStartupMessages(source(paste0(opt$codeDir, "CredibleSetTools.R")))

##############################################################################
## Load common data

# convert all boolean to caps strings
## TODO: Change these to booleans instead of string!
opt$cellType <- as.logical(toupper(opt$cellType))
opt$TargetGene <- as.logical(toupper(opt$TargetGene))
opt$TargetGeneTSS <- as.logical(toupper(opt$TargetGeneTSS))

# All genes
## TODO: Move all of this gene processing logic somewhere else, so that the gene files input into this script are fully ready
genes <- readBed(opt$genes)
genes$symbol <- unlist(lapply(strsplit(as.character(as.matrix(genes$name)), ";"), "[", 1))
genes <- addTSSToBED(genes)

# Genes used in the predictions
genes.uniq <- readBed(opt$genesUniq)
genes.uniq <- addTSSToBED(genes.uniq)

## Only include protein-coding genes
genes <- subset(genes, !grepl("NR_",name) & symbol %in% genes.uniq$name)
genes.uniq <- subset(genes.uniq, name %in% genes$symbol)
  
#  Vector of all diseases/traits in the variant list
variant.list <- read.delim(opt$variants, check.names=F)


# Finding a vector of relevant cell types
if (!is.null(opt$cellTypeTable)) {  
  cell.type.annot <- read.delim(opt$cellTypeTable, check.names=F)
  cell.type.list <- cell.type.annot$CellType
} else {
  if (opt$cellType) {
    ## JME: This is not the best way to get this list of cell types
    data <- read.delim(opt$bgOverlap, check.names=F, header=F)
    cell.type.list <- unique(data$V8)
  }
  # Assuming the 8th column is where the CellType is 
  else {
    ## JME: I don't understand what this does?
    cell.type.list <- "CELLTYPE"
  }
}


# All credible sets
all.cs <- read.delim(opt$credibleSets, check.names=F, stringsAsFactors=F)
all.cs$CredibleSet <- factor(all.cs$CredibleSet)
all.cs$MaxVariantScore <- sapply(all.cs$CredibleSet, function(cs) max(subset(variant.list, CredibleSet == cs)[,opt$variantScoreCol]))

# Optionally removing all coding, splice site, and promoter variants. These sets
# are used in the enrichment analysis.
filter.cs <- subset(all.cs, (!AnyCoding | !opt$removeCodingVariants) &
                      (!AnySpliceSite | !opt$removeSpliceSiteVariants) &
                      (!AnyPromoter | !opt$removePromoterVariants))  ## Consider whether to factor this out of this script

variant.list.filter <- subset(variant.list, CredibleSet %in% filter.cs$CredibleSet)

# Finding CSs and variants that exceed the provided threhold. If no variant score
# and threshold are provided, using all variants in downstream steps
sigScore.cs = NULL
variant.list.sigScore = NULL
if (!(is.null(opt$variantScoreCol)) & !(is.null(opt$variantScoreThreshold))) {
  sigScore.cs <- subset(filter.cs, CredibleSet %in% subset(variant.list.filter, get(opt$variantScoreCol) >= opt$variantScoreThreshold)$CredibleSet)
  variant.list.sigScore <- subset(variant.list.filter, CredibleSet %in% sigScore.cs$CredibleSet) ## this is the list of variants in credible sets with at least 1 variant with posterior prob >10%
}

# Overlapping variants with predictions
# the required column names
overlap <- loadVariantOverlap(opt$predictionFile, genes.uniq, genes, variant.names=variant.list$variant, isTargetGeneTSS=opt$TargetGeneTSS)
overlap <- filterVariantOverlap(overlap, opt$predScoreCol, opt$minPredScore, opt$minPredScorePromoters)


# Annotating overlaps
all.flat <- annotateVariantOverlaps(overlap, variant.list, all.cs)
if (opt$cellType) {
  all.flat <- subset(all.flat, CellType %in% cell.type.list)  ## IMPORTANT CHANGE
}
filter.flat <- subset(all.flat, CredibleSet %in% filter.cs$CredibleSet)

# If a threshold is provided, getting annotations for significant variants
sigScore.flat = NULL
if (!(is.null(opt$variantScoreCol)) & !(is.null(opt$variantScoreThreshold))) {
  sigScore.flat <- subset(filter.flat, CredibleSet %in% sigScore.cs$CredibleSet)
}
# Finding all diseases/traits in the credible set table
# Assume that the cs has only one trait
#traits <- unique(all.cs$Disease)
trait <- opt$trait

# Creating an output directory and writing the result to files
dir.create(paste0(opt$outbase,"data/"))
write.tab(all.flat, file=paste0(opt$outbase, "data/all.flat.tsv"))
write.tab(filter.flat, file=paste0(opt$outbase, "data/filter.flat.tsv"))
if (!(is.null(sigScore.flat))) {
  write.tab(sigScore.flat, file=paste0(opt$outbase, "data/sigScore.flat.tsv"))
}


##############################################################################
## Calculate enrichment per cell type by comparing the significant/provided
## variants with a set of background variants

# A set of background variants
bgVars <- read.delim(gzfile(opt$backgroundVariants), check.names=F, header=F)  ## TODO: Only need the count from this file
bgOverlap <-read.delim(gzfile(opt$bgOverlap), check.names=F, header=F)  ## TODO: The prediction - bgOverlap is the same for every trait


edir <- paste0(opt$outbase, "enrichment/")
dir.create(edir)


variant.by.cells <- getVariantByCellsTable(filter.flat, isTargetGene=opt$TargetGene, isCellType=opt$cellType)
#write.tab(variant.by.cells, file="variant.by.cells.tsv")

# With promoters
enrich <- computeCellTypeEnrichment(variant.by.cells,
                                    variant.list.filter,
                                    cell.type.list,
                                    trait,
                                    score.col=opt$variantScoreCol,
                                    min.score=opt$variantScoreThreshold,
                                    bg.vars=bgVars,
                                    bg.overlap=bgOverlap,
                                    isCellType=opt$cellType,
                                    enrichment.threshold=opt$biosampleEnrichThreshold)



print("Saving")
write.tab(enrich, file=opt$outEnrichment)


####################################################################################
## Write gene predictions table

enriched.cell.types <- subset(enrich, Significant)$CellType

gp.all <- getGenePrioritizationTable(
  all.flat, 
  all.cs, 
  genes, 
  genes.uniq, 
  enriched.cell.types, 
  cell.type.annot, 
  score.col=opt$predScoreCol, 
  score.min=-Inf,
  var.score.col=opt$variantScoreCol, 
  var.score.min=opt$variantScoreThreshold,
  max.distance=opt$genePredMaxDistance,
  method.name=opt$methodName)

writeGenePrioritizationTable(gp.all, file=paste0(opt$outGenePredTable))


