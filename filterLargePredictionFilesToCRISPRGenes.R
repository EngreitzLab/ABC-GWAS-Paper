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


#AllCellTypes
all.cell.expt <- fread('/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/comparisonRuns/AllCellTypes-ABC_comparison/experimentalData/experimentalData.AllCellTypes.txt')
cell.mapping <- as.data.table(list(expt.cell = c('GM12878','PrimaryHepatocytes','BJAB','BJAB_anti-IgM_anti-CD40_4hr','Jurkat','Jurkat_anti-CD3_PMA_4hr','THP1','THP1_LPS_4hr','NCCIT','LNCAP'),
                                   pred.cell = c('GM12878-Roadmap','liver-ENCODE','BJAB-Engreitz','BJAB_anti-IgM_anti-CD40_4hr-Engreitz','Jurkat-Engreitz','Jurkat_anti-CD3_PMA_4hr-Engreitz','THP1-Engreitz','THP1_LPS_4hr-Engreitz','NCCIT','LNCAP')))

lapply(cell.mapping$expt.cell, function (s) {
  pred.name <- subset(cell.mapping, expt.cell == s)$pred.cell
  filter_file(paste0("/seq/lincRNA/RAP/ABC/191216_ABC/ABC_out/",pred.name,"/Predictions_AvgHiC/EnhancerPredictionsAllPutative.bedpe.gz"), 
              s, 
              all.cell.expt, 
              paste0("/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/predictions/ABC.",pred.name,".AvgHiC.EnhancerPredictionsAllPutative.OnlyCRISPRGenes.bedpe"))
})
  

#Roadmap Links
# filter_file("/seq/lincRNA/RAP/External/EP_Predictions/Liu2017-Roadmap/Liu2017-EnhancerGeneLinksMerged.tsv.gz", 
#             "K562", 
#             k562.expt, 
#             "/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/predictions/Roadmap-Liu2017-EnhancerGeneLinksMerged.OnlyCRISPRGenes.tsv")
roadmap.links <- fread('/seq/lincRNA/RAP/External/EP_Predictions/Liu2017-Roadmap/Liu2017-EnhancerGeneLinksMerged.tsv.gz')
roadmap.links.neededForPaper <- subset(roadmap.links, (CellType == 'K562_Leukemia' & GeneSymbol %in% unique(subset(all.cell.expt, CellType == 'K562')$GeneSymbol)) | #K562
                                                      (CellType == 'Monocytes-CD14+_RO01746' & GeneSymbol %in% unique(subset(all.cell.expt, CellType == 'THP1')$GeneSymbol)) | #THP1
                                                      (CellType == 'CD19_Primary_Cells_Peripheral_UW' & GeneSymbol %in% unique(subset(all.cell.expt, CellType == 'BJAB')$GeneSymbol)) | #BJAB
                                         (CellType == 'CD3_Primary_Cells_Peripheral_UW' & GeneSymbol %in% unique(subset(all.cell.expt, CellType == 'Jurkat')$GeneSymbol)) | #Jurkat
                                         (CellType == 'GM12878_Lymphoblastoid' & GeneSymbol %in% unique(subset(all.cell.expt, CellType == 'GM12878')$GeneSymbol)) #GM
                                       )
write.table(roadmap.links.neededForPaper, "/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/predictions/Roadmap-Liu2017-EnhancerGeneLinksMerged.OnlyCRISPRGenes.tsv",
            sep = '\t', quote = F, col.names = T, row.names = F)

subset(roadmap.links, CellType == 'K562_Leukemia' & GeneSymbol == 'HBE1' & abs(startElement - 5497365) < 1000 )
subset(roadmap.links.neededForPaper, CellType == 'K562_Leukemia' & GeneSymbol == 'HBE1' & abs(startElement - 5497365) < 1000 )

#Roadmap
roadmap.expt <- fread('/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/comparisonRuns/RoadmapComparison-RoadmapMatched/experimentalData/experimentalData.RoadmapComparison-RoadmapMatched.txt')
cell.mapping <- as.data.table(list(expt.cell = c('BJAB','Jurkat','THP1'),
                                   pred.cell = c('CD19-positive_B_cell-Roadmap','CD3-positive_T_cell-Roadmap','CD14-positive_monocytes-Roadmap')))

lapply(cell.mapping$expt.cell, function (s) {
  pred.name <- subset(cell.mapping, expt.cell == s)$pred.cell
  filter_file(paste0("/seq/lincRNA/RAP/ABC/191216_ABC/ABC_out/",pred.name,"/Predictions_AvgHiC/EnhancerPredictionsAllPutative.bedpe.gz"), 
              s, 
              roadmap.expt, 
              paste0("/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/predictions/ABC.",pred.name,".AvgHiC.EnhancerPredictionsAllPutative.OnlyCRISPRGenes.bedpe"))
})

#PCHIC comparison
pchic.expt <- fread('/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/comparisonRuns/PCHiC-Comparison/experimentalData/experimentalData.PCHiC-Comparison.txt')
filter_file("//seq/lincRNA/RAP/ABC/191216_ABC/ABC_out/erythroblast-Corces2016//Predictions_AvgHiC/EnhancerPredictionsAllPutative.bedpe.gz", 
            "K562", 
            pchic.expt, 
            "/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/predictions/ABC/ABC.erythroblast-Corces2016.AvgHiC.EnhancerPredictionsAllPutative.OnlyCRISPRGenes.bedpe")

#PCHiC comparison - using PCHIC
cell.mapping <- as.data.table(list(expt.cell = c('BJAB','Jurkat','THP1','K562'),
                                   pred.cell = c('CD19-positive_B_cell-Roadmap','CD3-positive_T_cell-Roadmap','CD14-positive_monocytes-Roadmap','erythroblast-Corces2016')))

lapply(cell.mapping$expt.cell[[4]], function (s) {
  pred.name <- subset(cell.mapping, expt.cell == s)$pred.cell
  filter_file(paste0("/seq/lincRNA/RAP/ABC/191216_ABC/ABC_out/",pred.name,"/Predictions/EnhancerPredictionsAllPutative.bedpe.gz"), 
              s, 
              pchic.expt, 
              paste0("/seq/lincRNA/RAP/GWAS/210223_ABCGWASPaper_github/ABC-GWAS-Paper/comparePredictorsToCRISPRData/predictions/ABC/ABC.",pred.name,".PCHiC.EnhancerPredictionsAllPutative.OnlyCRISPRGenes.bedpe"))
})
