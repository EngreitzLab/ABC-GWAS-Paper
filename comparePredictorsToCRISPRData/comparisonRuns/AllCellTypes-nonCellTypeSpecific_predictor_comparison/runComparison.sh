Rscript ../../code/comparePredictionsToExperiment.R \
                --predictions config/pred.table.listing.txt \
                --experimentalData experimentalData/experimentalData-AllCellTypes_nonCellTypeSpecific.txt \
                --experimentalPositiveColumn 'Regulated' \
                --plotConfig config/plot.config.txt \
                --predConfig ../../config/pred.config.txt \
                --code ../../code/comparison.R \
                --outDir out/