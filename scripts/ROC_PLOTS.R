source("./scripts/checkLibraries.R")

plotROC_without_zigzag <- function(outputfolder, species_tissue, infoSamples, FDRInverse_cutoff, bh, bh_auc, fdr_inverse, fdr_inverse_auc){
 
  pdf(paste0(outputfolder, "/ROC_curve_",species_tissue,"_",infoSamples, ".pdf"), width = 10, height = 6)
  ### Extremum Distance Estimator for inflection point
  A = inflection::ede(as.numeric(bh$FPR), as.numeric(bh$TPR), 0)
  B = inflection::ede(as.numeric(fdr_inverse$FPR), as.numeric(fdr_inverse$TPR), 0)
  plot(bh$FPR, bh$TPR, col="red",  axes=F, main="ROC", type="l", lwd=3, ylab="TPR", xlab="FPR")
  mtext(species_tissue)
  axis(side=2)
  axis(side=1)
  lines(fdr_inverse$FPR, fdr_inverse$TPR, col="darkblue", lwd=3)
  lines(x = c(0,1), y = c(0,1))
  

  bh$valueTreshold <- as.numeric(bh$valueTreshold)
  bh$TPR <- as.numeric(bh$TPR)
  bh$FPR <- as.numeric(bh$FPR)
  bh <- na.omit(bh)
  
  fdr_inverse$valueTreshold <- as.numeric(fdr_inverse$valueTreshold)
  fdr_inverse$TPR <- as.numeric(fdr_inverse$TPR)
  fdr_inverse$FPR <- as.numeric(fdr_inverse$FPR)
  fdr_inverse <- na.omit(fdr_inverse)
  
  
  ## bh inflection point
  points(A[,3], min(bh$TPR[which(bh$FPR >= A[,3])]), lwd=5, pch=20, col="black")
  if(species_tissue == "Human lung"){
    points(B[,3], max(fdr_inverse$TPR[which(fdr_inverse$FPR <= B[,3])])+0.02, lwd=5, pch=20, col="black")
  } else{
    points(B[,3], max(fdr_inverse$TPR[which(fdr_inverse$FPR <= B[,3])]), lwd=5, pch=20, col="black")
  }
  
  ### add point based on threshold cut-off--> bh threshold cut-off <= 0.05
  points(max(bh$FPR[which(bh$valueTreshold <= 0.05)]), max(bh$TPR[which(bh$valueTreshold <= 0.05)]), lwd=2, pch=10, col="red")
  ## frd_inverse threshold cut-off <=  (that represents 5% of the FDR of all samples)
  points(max(fdr_inverse$FPR[which(fdr_inverse$valueTreshold <= FDRInverse_cutoff)]), max(fdr_inverse$TPR[which(fdr_inverse$valueTreshold <= FDRInverse_cutoff)]), lwd=2, pch=10, col="darkblue")
  legend("bottomright", legend = c(paste0("BH = ",  round(bh_auc,4)), paste0("FDR_inverse = ", round(fdr_inverse_auc,4)), " ", "inflection points", expression("pValue "<= 0.05), expression("qValue " <= 0.05)), col = c("red", "darkblue", "white","black", "red", "darkblue"), lty=c(1,1,NA,NA,NA,NA,NA),pch = c(NA,NA,NA,20,10,10,10), bty = "n", title = "AUC",inset = c(0.001, 0.35), lwd=3)
   dev.off()
   

   #### export stats  
   ## return thresholds (pValue or qValue) of inflection points
   pValue_inflection <- min(bh$valueTreshold[which(bh$FPR >= A[,3])])
   getAllRow_InflectionPointBH <- bh[bh$valueTreshold==pValue_inflection,]
   
   fdr_inverse_inflection <- max(fdr_inverse$valueTreshold[which(fdr_inverse$FPR <= B[,3])])
   getAllRow_InflectionPoint_FDR <- fdr_inverse[fdr_inverse$valueTreshold==fdr_inverse_inflection,]
   
   ## return thresholds (pValue or qValue) of thresholds
   pValue_threshold <- max(bh$valueTreshold[which(bh$valueTreshold <= 0.05)])
   getAllRow_ThresholdBH <- bh[bh$valueTreshold == pValue_threshold,]
   
   getAllRow_ThresholdFDR <- fdr_inverse[fdr_inverse$valueTreshold <= FDRInverse_cutoff,]
   getAllRow_ThresholdFDR <- fdr_inverse[fdr_inverse$valueTreshold == max(getAllRow_ThresholdFDR$valueTreshold),]
   
   collectAllStats <- list(pValue_inflection, getAllRow_InflectionPointBH, fdr_inverse_inflection, getAllRow_InflectionPoint_FDR, pValue_threshold, getAllRow_ThresholdBH, FDRInverse_cutoff, getAllRow_ThresholdFDR)
   names(collectAllStats) <- c("Inflection point BH", "TPR and FPR for inflection point BH", "Inflection point FDR_inverse", "TPR and FPR for inflection point FDR_inverse", 
                               "BH threshold of pValue < 0.05",  "TPR and FPR for BH threshold of pValue < 0.05", "FDR_inverse cut off for 5% of false discovery rate", "TPR and FPR for a cut off using 5% of false discovery rate" )
   capture.output(collectAllStats, file = file.path(outputfolder,paste0("Statistics_",species_tissue, "_",infoSamples,".tsv")))
 
}


plotROC_all_methods <- function(outputfolder, species_tissue, infoSamples, FDRInverse_cutoff, bh, bh_auc, fdr_inverse, fdr_inverse_auc, zigzag, zigzag_auc){
  
  pdf(paste0(outputfolder, "/ROC_curve_",species_tissue,"_",infoSamples, ".pdf"), width = 10, height = 6)
  
  ### Extremum Distance Estimator for inflection point
  A = inflection::ede(as.numeric(bh$FPR), as.numeric(bh$TPR), 0)
  B = inflection::ede(as.numeric(fdr_inverse$FPR), as.numeric(fdr_inverse$TPR), 0)
  C = inflection::ede(as.numeric(zigzag$FPR), as.numeric(zigzag$TPR), 0)
  plot(bh$FPR, bh$TPR, col="red",  axes=F, main="ROC", type="l", lwd=3, ylab="TPR", xlab="FPR")
  mtext(species_tissue)
  axis(side=2)
  axis(side=1)
  lines(fdr_inverse$FPR, fdr_inverse$TPR, col="darkblue", lwd=3)
  lines(zigzag$FPR, zigzag$TPR, col="gray", lwd=3)
  lines(x = c(0,1), y = c(0,1))
 
  bh$valueTreshold <- as.numeric(bh$valueTreshold)
  bh$TPR <- as.numeric(bh$TPR)
  bh$FPR <- as.numeric(bh$FPR)
  bh <- na.omit(bh)
  
  fdr_inverse$valueTreshold <- as.numeric(fdr_inverse$valueTreshold)
  fdr_inverse$TPR <- as.numeric(fdr_inverse$TPR)
  fdr_inverse$FPR <- as.numeric(fdr_inverse$FPR)
  fdr_inverse <- na.omit(fdr_inverse)
  
  zigzag$valueTreshold <- as.numeric(zigzag$valueTreshold)
  zigzag$TPR <- as.numeric(zigzag$TPR)
  zigzag$FPR <- as.numeric(zigzag$FPR)
  zigzag <- na.omit(zigzag)
  
  ## bh inflection point
  points(A[,3], min(bh$TPR[which(bh$FPR >= A[,3])]), lwd=5, pch=20, col="black")
  if(species_tissue == "Human lung"){
    points(B[,3], max(fdr_inverse$TPR[which(fdr_inverse$FPR <= B[,3])])+0.02, lwd=5, pch=20, col="black")
  } else{
    points(B[,3], max(fdr_inverse$TPR[which(fdr_inverse$FPR <= B[,3])]), lwd=5, pch=20, col="black")
  }
  if(species_tissue == "Drosophila testis"){
    points(C[,3]+0.008, min(zigzag$TPR[which(zigzag$FPR >= C[,3])]), lwd=5, pch=20, col="black")
  } else {
    ## zigzag inflection point
    points(C[,3], min(zigzag$TPR[which(zigzag$FPR >= C[,3])]), lwd=5, pch=20, col="black")
  }
  
  ### add point based on threshold cut-off--> bh threshold cut-off <= 0.05
  points(max(bh$FPR[which(bh$valueTreshold <= 0.05)]), max(bh$TPR[which(bh$valueTreshold <= 0.05)]), lwd=2, pch=10, col="red")
  ## frd_inverse threshold cut-off <=  (that represents 5% of the FDR of all samples)
  points(max(fdr_inverse$FPR[which(fdr_inverse$valueTreshold <= FDRInverse_cutoff)]), max(fdr_inverse$TPR[which(fdr_inverse$valueTreshold <= FDRInverse_cutoff)]), lwd=2, pch=10, col="darkblue")
  ## zigzag threshold cut-off >= 0.95
  points(max(zigzag$FPR[which(zigzag$valueTreshold >= 0.95)]), max(zigzag$TPR[which(zigzag$valueTreshold >= 0.95)]), lwd=2, pch=10, col="gray")
  legend("bottomright", legend = c(paste0("BH = ",  round(bh_auc,4)), paste0("FDR_inverse = ", round(fdr_inverse_auc,4)), paste0("zigzag = ", round(zigzag_auc,4)), " ", "inflection points", expression("pValue "<= 0.05), expression("qValue " <= 0.05), expression("PP ">= 0.95)), col = c("red", "darkblue", "gray", "white","black", "red", "darkblue", "gray"), lty=c(1,1,1,NA,NA,NA,NA,NA),pch = c(NA,NA,NA,NA,20,10,10,10), bty = "n", title = "AUC",inset = c(0.001, 0.35), lwd=3)
  dev.off()
  
  
  #### export stats  
  ## return thresholds (pValue or qValue) of inflection points
  pValue_inflection <- min(bh$valueTreshold[which(bh$FPR >= A[,3])])
  getAllRow_InflectionPointBH <- bh[bh$valueTreshold==pValue_inflection,]
  
  fdr_inverse_inflection <- max(fdr_inverse$valueTreshold[which(fdr_inverse$FPR <= B[,3])])
  getAllRow_InflectionPoint_FDR <- fdr_inverse[fdr_inverse$valueTreshold==fdr_inverse_inflection,]
  
  getAllRow_InflectionPoint_zigzag <- zigzag[zigzag$TPR==min(zigzag$TPR[which(zigzag$FPR >= C[,3])]),]
  

  
  ## return thresholds (pValue or qValue) of thresholds
  pValue_threshold <- max(bh$valueTreshold[which(bh$valueTreshold <= 0.05)])
  getAllRow_ThresholdBH <- bh[bh$valueTreshold == pValue_threshold,]
  
  getAllRow_ThresholdFDR <- fdr_inverse[fdr_inverse$valueTreshold <= FDRInverse_cutoff,]
  getAllRow_ThresholdFDR <- fdr_inverse[fdr_inverse$valueTreshold == max(getAllRow_ThresholdFDR$valueTreshold),]
  
  getAllRow_ThresholdZigzag <- zigzag[zigzag$TPR==max(zigzag$TPR[which(zigzag$valueTreshold >= 0.95)]),]
  
  
  collectAllStats <- list(pValue_inflection, getAllRow_InflectionPointBH, fdr_inverse_inflection, getAllRow_InflectionPoint_FDR, getAllRow_InflectionPoint_zigzag,
                          pValue_threshold, getAllRow_ThresholdBH, FDRInverse_cutoff, getAllRow_ThresholdFDR, getAllRow_ThresholdZigzag)
  names(collectAllStats) <- c("Inflection point BH", "TPR and FPR for inflection point BH", "Inflection point FDR_inverse", "TPR and FPR for inflection point FDR_inverse", "TPR and FPR for inflection point zigzag",
                              "BH threshold of pValue < 0.05",  "TPR and FPR for BH threshold of pValue < 0.05", "FDR_inverse cut off for 5% of false discovery rate", "TPR and FPR for a cut off using 5% of false discovery rate",
                              "TPR and FPR for a threshold >= 0.95")
  capture.output(collectAllStats, file = file.path(outputfolder,paste0("Statistics_",species_tissue, "_",infoSamples,".tsv")))
  
  
 
}