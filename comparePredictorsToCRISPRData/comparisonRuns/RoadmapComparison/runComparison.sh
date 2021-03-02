Rscript ../../code/comparePredictionsToExperiment.R \
                --predictions config/pred.table.listing.txt \
                --experimentalData experimentalData/experimentalData.RoadmapComparison.txt \
                --experimentalPositiveColumn 'Regulated' \
                --cellNameMapping config/cellNameMapping.txt \
                --plotConfig config/plot.config.txt \
                --predConfig ../../config/pred.config.txt \
                --code ../../code/comparison.R \
                --outDir out/