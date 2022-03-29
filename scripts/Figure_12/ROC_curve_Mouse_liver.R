source("./scripts/ROC_AUC.R")
source("./scripts/ROC_PLOTS.R")

TP_TN_file <- "./data/benchmark_files/TP_TN_usingJust_ProteinCoding_TP_filter_TN_others.tsv"

bh_file_89Samples <- "./data/Combining_process/Mouse/Calls_merging_BH_cutoff=0.05_species_id=10090_89Libraries.tsv"
fdr_inverse_file_89Samples <- "./data/Combining_process/Mouse/Calls_merging_fdr_inverse_cutoff=0.05_species_id=10090_89Libraries.tsv"
zigzag_file_89Samples <- "./data/Combining_process/Mouse/Calls_zigzag_89Libraries.tsv"

outputFolder <- "./figures/Figure_12/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}


## Plot ROC using all 89 samples for all methods
collectData <- collectInfoMethods(TP_TN_file=TP_TN_file, bh_file=bh_file_89Samples, fdr_inverse_file=fdr_inverse_file_89Samples, zigzag_file=zigzag_file_89Samples)

selectTable <- collectData[[2]]
selectTable <- na.omit(selectTable) 
bh <- ROC(table = selectTable, method = "bh")
bh <- bh[!grepl("NA", bh$valueTreshold),]
bh_auc <- trapezoid(x=as.numeric(bh$FPR), as.numeric(bh$TPR))

selectTable <- collectData[[3]]
fdr_inverse <- ROC(table = selectTable, method = "fdr_inverse")
fdr_inverse <- fdr_inverse[order(as.numeric(fdr_inverse$TPR)),]
fdr_inverse_auc <- trapezoid(x=as.numeric(fdr_inverse$FPR), as.numeric(fdr_inverse$TPR))

selectTable <- collectData[[4]]
zigzag <- ROC(table = selectTable, method = "zigzag")
zigzag <- zigzag[!grepl("NA", zigzag$valueTreshold),]
zigzag <- zigzag[order(as.numeric(zigzag$TPR)),]
zigzag_auc <- trapezoid(x=as.numeric(zigzag$FPR), as.numeric(zigzag$TPR))

plotROC_all_methods(outputfolder = outputFolder, species_tissue = "Mouse liver",  infoSamples = "89_Samples", FDRInverse_cutoff = collectData[[1]], bh = bh, bh_auc = bh_auc, fdr_inverse = fdr_inverse, fdr_inverse_auc = fdr_inverse_auc, zigzag = zigzag, zigzag_auc = zigzag_auc)
