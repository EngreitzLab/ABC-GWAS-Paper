### Code to create PR curves (Extended Data Fig 4)

Each panel in Extended Data Fig. 4 represents a comparison of multiple predictive models against a relevant subset of the experimental CRISPR dataset. This is a self contained repository which contains the relevant (subsets of) the necessary prediction files and the experimental data files.

PR curves for each comparison are available at:`comparisonRuns/*/out/combined/pr.curve.txt`. The results in the `comparisonRuns/*/out/` directories can be produced by running the script `comparisonRuns/*/runComparison.sh`. See [here](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/tree/better_comparison/comparison) for more documentation on the code which produces the comparison.

Each `comparisonRuns/` directory contains the appropriate experimental data file for this comparison (located in `comparisonRuns/*/experimentalData/`. This file is equivalent to the first 18 columns of Supplementary Table 5 which is further subsetted to the appropriate cell types for the comparison. It also contains pointers to the necessary prediction files (located at `comparisonRuns/*/config/pred.table.listing.txt`). The ABC prediction files located in `predictions/ABC/` contain the ABC predictions needed for these PR curves. These files only contain predictions for genes for which CRISPR data is available. The full set of ABC predictions is available on the FTP site: `http://ftp.broadinstitute.org/outgoing/lincRNA/ABC/Full-ABC-Output/`

Some figures in the manuscript contain PR curves from multiple comparisonRuns folders - the relationship is as follows:

* ED 4a: `K562-only`
* ED 4c: `AllCellTypes-ABC_comparison` (contains red ABC pr curve), `AllCellTypes-nonCellTypeSpecific_predictor_comparison` (contains all other pr curves)
* ED 4d: `RoadmapComparison-RoadmapMatched` (contains green PR curve), `RoadmapComparison` (contains other pr curves)
* ED 4e: `PCHiC-Comparison`


