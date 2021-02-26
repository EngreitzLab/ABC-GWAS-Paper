library(data.table)

#Subset ABC predictions files to only include genes we have CRISPR data for. 
#This is to make the make the predictions files small enough that they can be placed on github


filter_file <- function(pred_file, expt_cell_type, expt, new_file) {
  pred <- fread(pred_file)
  genes <- unique(subset(expt, CellType == expt_cell_type)$GeneSymbol)
  pred <- subset(pred, GeneSymbol %in% genes)

  write.table(pred, new_file, sep = '\t', quote = F, col.names = T, row.names = F)
}

#K562-only
k562.expt <- fread("/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/comparisonRuns/K562-only/experimentalData/experimentalData.K562-only.txt")
filter_file("//seq/lincRNA/RAP/ABC/191216_ABC/ABC_out/K562-Roadmap/Predictions_AvgHiC/EnhancerPredictionsAllPutative.bedpe.gz", 
            "K562", 
            k562.expt, 
            "/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/predictions/ABC.K562.AvgHiC.EnhancerPredictionsAllPutative.OnlyCRISPRGenes.bedpe")

filter_file("/seq/lincRNA/RAP/External/EP_Predictions/Liu2017-Roadmap/Liu2017-EnhancerGeneLinksMerged.tsv.gz", 
            "K562", 
            k562.expt, 
            "/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/predictions/Roadmap-Liu2017-EnhancerGeneLinksMerged.OnlyCRISPRGenes.tsv")

#AllCellTypes
all.cell.expt <- fread('/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/comparisonRuns/AllCellTypes-ABC_comparison/experimentalData/experimentalData.AllCellTypes.txt')
cell.mapping <- as.data.table(list(expt.cell = c('GM12878','PrimaryHepatocytes','BJAB','BJAB_anti-IgM_anti-CD40_4hr','Jurkat','Jurkat_anti-CD3_PMA_4hr','THP1','THP1_LPS_4hr','NCCIT','LNCAP'),
                                   pred.cell = c('GM12878-Roadmap','liver-ENCODE','BJAB-Engreitz','BJAB_anti-IgM_anti-CD40_4hr-Engreitz','Jurkat-Engreitz','Jurkat_anti-CD3_PMA_4hr-Engreitz','THP1-Engreitz','THP1_LPS_4hr-Engreitz','NCCIT','LNCAP')))

lapply(cell.mapping$expt.cell[9:10], function (s) {
  pred.name <- subset(cell.mapping, expt.cell == s)$pred.cell
  filter_file(paste0("/seq/lincRNA/RAP/ABC/191216_ABC/ABC_out/",pred.name,"/Predictions_AvgHiC/EnhancerPredictionsAllPutative.bedpe.gz"), 
              s, 
              all.cell.expt, 
              paste0("/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/predictions/ABC.",pred.name,".AvgHiC.EnhancerPredictionsAllPutative.OnlyCRISPRGenes.bedpe"))
})
  
