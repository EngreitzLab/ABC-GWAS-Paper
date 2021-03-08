# coding: utf-8
from os.path import join
import pandas as pd


#configfile: "ABC-Max.config.json"  ## Specify this on the command line

pred_config = config["predictionsTable"]
trait_config = config["traitTable"]

preds_config_file = pd.read_table(pred_config).set_index("entry", drop=False)
trait_config_file = pd.read_table(trait_config).set_index("entry", drop=False)


rule all:
	input:
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapAllSNPs.tsv.gz"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.tsv"), pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.AllNoncoding.tsv"), pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/{trait}.variants.bed", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/{trait}.variants.bedgraph", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/{trait}.{pred}.tsv.gz", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{trait}/{trait}.{pred}.txt"), trait=config["traits"], pred=config["predictions"]),
		expand(os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.pdf"), trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/GenePredictions.allCredibleSets.tsv", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"]),
		expand("{outdir}{pred}/{trait}/GenePrecisionRecall.pdf", outdir=config["outDir"], trait=config["traits"], pred=config["predictions"])


rule computeBackgroundOverlap:
	input:
		predFile = lambda wildcard: preds_config_file.loc[wildcard.pred, "predFile"],
		allVariants = config["bgVariants"],
		chrSizes = config["chrSizes"],
		CDS = config["CDS"]
	params:
		cellType = lambda wildcard: str(preds_config_file.loc[wildcard.pred, "cellType"]), 
		outDir = expand("{outdir}{{pred}}", outdir=config["outDir"])
	output:
		overallOverlap = os.path.join(config["outDir"], "{pred}/{pred}.OverlapAllSNPs.tsv.gz"),
		overallOverlapCounts = os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.tsv"),
		noncodingOverlap = os.path.join(config["outDir"], "{pred}/{pred}.OverlapCounts.AllNoncoding.tsv")
	priority: 5
	log: os.path.join(config["logDir"], "{pred}.bgoverlap.log")
	message: "Overlapping background variants with predictions: {wildcards.pred}"
	run:
		shell(
			"""
			# TODO: find an alternative to deal with pipefail
			set +o pipefail;
			# Intersecting a background list of variants with predicted enhancers 
			# to compute background rate at which common variants overlap enhancers
			# make output dir
			if [ ! -d {params.outDir} ]
			then
				mkdir {params.outDir}
			fi
			# Compute fraction of variants overlapping predictions in each cell type
			# Finding the relevant columns
			if ({params.cellType}=='True')
			then
				zcat {input.predFile} | csvtk cut -t -f chr,start,end,CellType | sed 1d | sort -k 1,1 -k 2,2n | uniq | bedtools sort -i stdin -faidx {input.chrSizes} | \
				bedtools intersect -sorted -g {input.chrSizes} -a {input.allVariants} -b stdin -wa -wb | gzip > {output.overallOverlap};
			else
				zcat {input.predFile} | csvtk cut -t -f chr,start,end | sed 1d | sort -k 1,1 -k 2,2n | uniq | bedtools sort -i stdin -faidx {input.chrSizes} | \
				bedtools intersect -sorted -g {input.chrSizes} -a {input.allVariants} -b stdin -wa -wb | gzip > {output.overallOverlap};
			fi
			
			# Getting the cell type column and counting
			 zcat {output.overallOverlap} | cut -f 7 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\\t' > {output.overallOverlapCounts};

			# Compute fraction of noncoding variants overlapping predictions in any cell type
			 zcat {output.overallOverlap} | bedtools intersect -v -a stdin -b {input.CDS} | cut -f 1-3,7 | sort | uniq | cut -f 4 | sort | uniq -c | sed 's/^ *//' | tr ' ' '\\t' > {output.noncodingOverlap}			
			""")

rule createVarFiles:
	input:
		varList = lambda wildcard: trait_config_file.loc[wildcard.trait, "varList"]
	output:
		varBed = os.path.join(config["outDir"],"{pred}/{trait}/{trait}.variants.bed"),
		varBedgraph = os.path.join(config["outDir"],"{pred}/{trait}/{trait}.variants.bedgraph"),
		sigvarList = os.path.join(config["outDir"],"{pred}/{trait}/{trait}.sig.varList.tsv")
	log: os.path.join(config["logDir"], "{trait}.{pred}.createbed.log")
	priority: 4
	params:
		varFilterCol = lambda wildcard: trait_config_file.loc[wildcard.trait, "varFilterCol"],
		varFilterThreshold = lambda wildcard: trait_config_file.loc[wildcard.trait, "varFilterThreshold"],
		chrSizes = config["chrSizes"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/")
	message: "Creating variant BED files"
	run:
		if {params.varFilterCol} is not None:
			shell(
				"""
				# make output dir 
				if [ ! -d {params.outDir} ]
                       		then
                                	mkdir {params.outDir}
                        	fi
				# Subsetting the variant list based on significance
				# Finding the score colum
				#scoreCol=$(awk -v RS='\\t' '/{params.varFilterCol}/{{print NR; exit}}' {input.varList});

				# Filtering to retain variants exceeding the threshold
				#awk '{{ if ($($scoreCol) >= {params.varFilterThreshold}) {{ print }} }}' {input.varList}  > {output.sigvarList};
				cat {input.varList} | csvtk -t filter -f "{params.varFilterCol}>={params.varFilterThreshold}" > {output.sigvarList};
	
				# Creating the bed file
				# Finding and cutting chr, position, and variant columns
				# TODO: do not require start and stop, only position?
				cat {output.sigvarList} | csvtk cut -t -f chr,position,variant | sed '1d' | awk -F "\\t" "\$1 = \$1 FS \$2-1 FS \$2 FS \$3 FS" | cut -f1-4 | sed -e 's/8.1e+07/81000000/g'> {output.varBed};

				# Ensure that variants are sorted for bedtools -sorted overlap algorithm
				cat {output.varBed} | bedtools sort -i stdin -faidx {params.chrSizes} | uniq > {output.varBedgraph};
				""")
		else:
			shell(
				"""
				#fi
				# Creating the bed file for all variants
				# Finding and cutting chr, pos, and var columns
				cat {input.varList} | csvtk cut -t -f chr,position,variant | sed '1d' | awk -F "\\t" "\$1 = \$1 FS \$2-1 FS \$2 FS \$3 FS" | cut -f1-4 > {output.varBed};

				# Ensure that variants are sorted for bedtools -sorted overlap algorithm
				cat {output.varBed} | bedtools sort -i stdin -faidx {params.chrSizes} | uniq > {output.varBedgraph};
				""")
	

rule overlapVariants:
	input:
		predFile = lambda wildcard: preds_config_file.loc[wildcard.pred, "predFile"],
		varList = lambda wildcard: trait_config_file.loc[wildcard.trait, "varList"],
		varBed = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.variants.bed", outdir=config["outDir"]),
		varBedgraph = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.variants.bedgraph", outdir=config["outDir"])
	output:
		overlap = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.{{pred}}.tsv.gz", outdir=config["outDir"])
	priority: 3
	log: os.path.join(config["logDir"], "{trait}.{pred}.overlap.log")
	params:
		chrSizes = config["chrSizes"]
	message: "Overlapping {wildcards.trait} variants with {wildcards.pred} enhancers"
	run:
		shell(
			"""
			# TODO: find an alternative to deal with pipefail
			set +o pipefail;

			# Creating an empty file with the final columns
			zcat {input.predFile} | head -1 | awk '{{ print $0 "\\tvariant.chr\\tvariant.start\\tvariant.end\\tQueryRegionName" }}' | gzip > {output.overlap};
	
			# Intersecting variants with predictions
			zcat {input.predFile} | sed 1d | bedtools intersect -g {params.chrSizes} -b {input.varBedgraph} -a stdin -wb | gzip >> {output.overlap}
			""")


rule annotateVariants:
	input:
		predFile = lambda wildcard: preds_config_file.loc[wildcard.pred, "predFile"],
		varList = lambda wildcard: trait_config_file.loc[wildcard.trait, "varList"],
		csList = lambda wildcard: trait_config_file.loc[wildcard.trait, "csList"],
		predOverlapFile = expand("{outdir}{{pred}}/{{trait}}/{{trait}}.{{pred}}.tsv.gz", outdir=config["outDir"]),
		bgVars = config["bgVariants"],
		bgOverlap = expand("{outdir}{{pred}}/{{pred}}.OverlapAllSNPs.tsv.gz", outdir=config["outDir"]) 
	output:
		touch(os.path.join(config["outDir"], "{pred}/{trait}/{trait}.{pred}.txt")),
		enrichFile = expand("{outdir}{{pred}}/{{trait}}/enrichment/Enrichment.CellType.vsScore.{{trait}}.tsv", outdir=config["outDir"]),
		genePredTable = expand("{outdir}{{pred}}/{{trait}}/GenePredictions.allCredibleSets.tsv", outdir=config["outDir"])
	priority: 2
	log: os.path.join(config["logDir"], "{trait}.{pred}.annotate.log")
	params:
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred,"celltypeAnnotation"],
		codeDir = config["codeDir"],
		projectDir = config["projectDir"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
		predScoreCol = lambda wildcard: preds_config_file.loc[wildcard.pred,"predScoreCol"],
		minPredScore = lambda wildcard: preds_config_file.loc[wildcard.pred,"minPredScore"],
		minPredScorePromoter = lambda wildcard: preds_config_file.loc[wildcard.pred,"minPredScorePromoter"],
		varScoreCol = lambda wildcard: trait_config_file.loc[wildcard.trait,"varFilterCol"],
		varScoreType = lambda wildcard: trait_config_file.loc[wildcard.trait,"varScoreType"],
		varScoreThreshold = lambda wildcard: trait_config_file.loc[wildcard.trait,"varFilterThreshold"],
		genes = lambda wildcard:preds_config_file.loc[wildcard.pred,"genes"],
		genesUniq = lambda wildcard: preds_config_file.loc[wildcard.pred,"genesUniq"],
		cellType = lambda wildcard: str(preds_config_file.loc[wildcard.pred,"cellType"]),
		isTargetGene = lambda wildcard: str(preds_config_file.loc[wildcard.pred,"TargetGeneTSS"])
	message: "Annotating {wildcards.trait} variants with {wildcards.pred} predictions"
	run:
		shell(
                """
                Rscript {params.projectDir}/src/AnnotateCredibleSets.R \
                --variants {input.varList} \
                --credibleSets {input.csList} \
                --predictionFile {input.predOverlapFile} \
                --methodName {wildcards.pred} \
                --outbase {params.outDir} \
                --outEnrichment {output.enrichFile} \
                --outGenePredTable {output.genePredTable} \
                --predScoreCol {params.predScoreCol} \
                --minPredScore {params.minPredScore} \
                --minPredScorePromoters {params.minPredScorePromoter} \
                --backgroundVariants {input.bgVars} \
                --bgOverlap {input.bgOverlap} \
                --trait {wildcards.trait} \
                --codeDir {params.codeDir} \
                --variantScoreCol {params.varScoreCol} \
                --variantScoreThreshold {params.varScoreThreshold} \
                --cellTypeTable {params.cellTypeTable} \
                --genes {params.genes} \
                --genesUniq {params.genesUniq} \
                --cellType {params.cellType} \
                --TargetGeneTSS {params.isTargetGene}
				""")	


rule plotTraitEnrichment:
	input: 
		cellTypeEnrichments = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.tsv")
	output:
		outpdf = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.pdf"),
		outeps = os.path.join(config["outDir"], "{pred}/{trait}/CellTypeEnrichment.{trait}.eps")
	params:
		cellTypeTable = lambda wildcard: preds_config_file.loc[wildcard.pred, "celltypeAnnotation"],
		projectDir = config["projectDir"],
		outDir = os.path.join(config["outDir"], "{pred}/{trait}/"),
		#cellTypeEnrichments_noPromoter = os.path.join(config["outDir"], "{pred}/{trait}/enrichment/Enrichment.CellType.vsScore.{trait}.noPromoter.tsv"),
	 	isCellType = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"cellType"]), 
		#hasPromoterColumn = lambda wildcard: bool(preds_config_file.loc[wildcard.pred,"hasPromoter"])
	priority: 1
	message: "Running enrichment plots"
	run:
		shell(
			"""
			Rscript {params.projectDir}src/PlotCellTypeEnrichment.R \
			--outdir {params.outDir} \
			--outPdf {output.outpdf} \
			--outEps {output.outeps} \
			--cellTypes {params.cellTypeTable} \
			--cellTypeEnrichments {input.cellTypeEnrichments} \
			--codeDir {params.projectDir} \
			--trait {wildcards.trait} 
			""")


rule plotGenePrecisionRecall:
	input:
		genePredTable = expand("{outdir}{{pred}}/{{trait}}/GenePredictions.allCredibleSets.tsv", outdir=config["outDir"]),
		knownGenes = lambda wildcard: trait_config_file.loc[wildcard.trait, "knownGenes"]
	output:
		prPdf = os.path.join(config["outDir"], "{pred}/{trait}/GenePrecisionRecall.pdf")
	params:
		codeDir = config["codeDir"],
		projectDir = config["projectDir"]
	message: "Running precision-recall plot for {wildcards.trait} and {wildcards.pred} predictions"
	run:
		shell(
			"""
			Rscript {params.projectDir}src/PlotGenePrecisionRecall.R \
			--outPdf {output.prPdf} \
			--genePredTable {input.genePredTable} \
			--knownGenes {input.knownGenes} \
			--codeDir {params.projectDir}
			""")

