source("./scripts/ROC_AUC.R")

TP_TN_file <- "./data/benchmark_files/markers_fromZigzag_Drosophila.tsv"
bh_file_0.001 <- "./complementary_analysis/results_merging_BH_frdInverse/Drosophila/Calls_merging_BH_cutoff=0.001_species_id=7227.tsv"
fdr_inverse_file_0.001 <- "./complementary_analysis/results_merging_BH_frdInverse/Drosophila/Calls_merging_fdr_inverse_cutoff=0.001_species_id=7227.tsv"

bh_file_0.0001 <- "./complementary_analysis/results_merging_BH_frdInverse/Drosophila/Calls_merging_BH_cutoff=1e-04_species_id=7227.tsv"
fdr_inverse_file_0.0001 <- "./complementary_analysis/results_merging_BH_frdInverse/Drosophila/Calls_merging_fdr_inverse_cutoff=1e-04_species_id=7227.tsv"

outputFolder <- "./stats_info/merging_suppl_Info/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}


collectStats_merging <- function(outputFolder, bh, fdr_inverse, cutoff){
  
  bh$valueTreshold <- as.numeric(bh$valueTreshold)
  bh$TPR <- as.numeric(bh$TPR)
  bh$FPR <- as.numeric(bh$FPR)
  bh <- na.omit(bh)
  
  fdr_inverse$valueTreshold <- as.numeric(fdr_inverse$valueTreshold)
  fdr_inverse$TPR <- as.numeric(fdr_inverse$TPR)
  fdr_inverse$FPR <- as.numeric(fdr_inverse$FPR)
  fdr_inverse <- na.omit(fdr_inverse)
  
  ## return thresholds (pValue or qValue) of thresholds
  pValue_threshold <- max(bh$valueTreshold[which(bh$valueTreshold <= cutoff)])
  getAllRow_ThresholdBH <- bh[bh$valueTreshold == pValue_threshold,]
  
  getAllRow_ThresholdFDR <- fdr_inverse[fdr_inverse$valueTreshold <= collectData[[1]],]
  getAllRow_ThresholdFDR <- fdr_inverse[fdr_inverse$valueTreshold == max(getAllRow_ThresholdFDR$valueTreshold),]
  
  
  collectAllStats <- list(getAllRow_ThresholdBH, getAllRow_ThresholdFDR)
  names(collectAllStats) <- c(paste0("TPR and FPR for BH threshold of pValue < ", cutoff), paste0("TPR and FPR for a cut off using ", cutoff*100, "% of false discovery rate"))
  capture.output(collectAllStats, file = file.path(outputFolder,paste0("Statistics_Drosophila_testis_",cutoff,"_threshold_for_RefInt_methods.tsv")))

}

collectData <- collectInfoMethods(TP_TN_file=TP_TN_file, bh_file=bh_file_0.001, fdr_inverse_file=fdr_inverse_file_0.001, zigzag_file=NULL)

selectTable <- collectData[[2]]
selectTable <- na.omit(selectTable) 
bh <- ROC(table = selectTable, method = "bh")
bh <- bh[!grepl("NA", bh$valueTreshold),]
bh_auc <- trapezoid(x=as.numeric(bh$FPR), as.numeric(bh$TPR))

selectTable <- collectData[[3]]
fdr_inverse <- ROC(table = selectTable, method = "fdr_inverse")
fdr_inverse <- fdr_inverse[order(as.numeric(fdr_inverse$TPR)),]
fdr_inverse_auc <- trapezoid(x=as.numeric(fdr_inverse$FPR), as.numeric(fdr_inverse$TPR))

collectStats_merging(outputFolder = outputFolder, bh=bh, fdr_inverse=fdr_inverse, cutoff=0.001)


collectData <- collectInfoMethods(TP_TN_file=TP_TN_file, bh_file=bh_file_0.0001, fdr_inverse_file=fdr_inverse_file_0.0001, zigzag_file=NULL)

selectTable <- collectData[[2]]
selectTable <- na.omit(selectTable) 
bh <- ROC(table = selectTable, method = "bh")
bh <- bh[!grepl("NA", bh$valueTreshold),]
bh_auc <- trapezoid(x=as.numeric(bh$FPR), as.numeric(bh$TPR))

selectTable <- collectData[[3]]
fdr_inverse <- ROC(table = selectTable, method = "fdr_inverse")
fdr_inverse <- fdr_inverse[order(as.numeric(fdr_inverse$TPR)),]
fdr_inverse_auc <- trapezoid(x=as.numeric(fdr_inverse$FPR), as.numeric(fdr_inverse$TPR))

collectStats_merging(outputFolder = outputFolder, bh=bh, fdr_inverse=fdr_inverse, cutoff=0.0001)
