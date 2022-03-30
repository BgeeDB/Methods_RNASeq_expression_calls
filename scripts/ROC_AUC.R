source("./scripts/checkLibraries.R")

## ROC function
ROC <- function(table, method){
  
  collectInfo <- c()
  for (i in unique(table$active)) {
    print(i)
    valueTreshold <- i
    
    if(method == "bh" | method == "fdr_inverse"){
      selectTable$callApproach <- ifelse(table$active <= i, "present", "absent")
    } else {
      selectTable$callApproach <- ifelse(table$active >= i, "present", "absent")
    }
    
    ## present in control markers file and also in the method
    TP <- dplyr::filter(selectTable, callMarkers == "present" & callApproach == "present")
    TP <- nrow(TP)
    ## absent in control markers file and present in method
    FP <- dplyr::filter(selectTable, callMarkers == "absent" & callApproach == "present")
    FP <- nrow(FP)
    ## absent in the control markers file and absent in the method
    TN <- dplyr::filter(selectTable, callMarkers == "absent" & callApproach == "absent")
    TN <- nrow(TN)
    ## present in the control markers file and absent in the method
    FN <- dplyr::filter(selectTable, callMarkers == "present" & callApproach == "absent")
    FN <- nrow(FN)
    ## Calculate ratios
    TPR <- TP/(TP + FN)
    FPR <- FP/(FP+TN)
    TP_FP_Info <- data.frame(valueTreshold,TP,FP,TN,FN,TPR,FPR)
    collectInfo <- rbind(collectInfo,TP_FP_Info)
  }
  ## force ROC to start (0,0) and finish (1,1)
  top <- c("NA","NA","NA","NA","NA", 0,0)
  bottom <- c("NA","NA","NA","NA","NA", 1,1)
  collectInfo <- rbind(collectInfo, bottom, top)
  ## sort result by ordering FPR
  collectInfo <- collectInfo[order(as.numeric(collectInfo$FPR)),]
  return(collectInfo)
}

### AUC function
trapezoid <- function(x, y) sum(diff(x)*(y[-1]+y[-length(y)]))/2

### Collect input data function
collectInfoMethods <- function(TP_TN_file, bh_file, fdr_inverse_file, zigzag_file){
  TP_TN <- read.table(TP_TN_file, header=TRUE, sep="\t")
  colnames(TP_TN) <- c("id","callMarkers")
  
  if(is.null(zigzag_file) == TRUE) {
    bh <- read.table(bh_file, header=TRUE, sep="\t")
    bh <- bh[,1:2]
    colnames(bh) <- c("id", "active")
    
    fdr_inverse <- read.table(fdr_inverse_file, header=TRUE, sep="\t")
    ## collect value of qValue cutoff where samples were called present for qValue = 0.05 per individual sample
    FDRInverse_cutoff <- max(fdr_inverse$minimum_qValue[which(fdr_inverse$call == "present")])
    fdr_inverse <- fdr_inverse[,1:2]
    colnames(fdr_inverse) <- c("id", "active")
    
    selectTableBH <- merge(TP_TN, bh, by ="id")
    selectTableBH <- selectTableBH[order(selectTableBH$active) , ]
    
    selectTableFDR <- merge(TP_TN, fdr_inverse, by ="id")
    selectTableFDR <- selectTableFDR[order(selectTableFDR$active) , ]
    
    return(list(FDRInverse_cutoff, selectTableBH, selectTableFDR))
    
  } else {
    bh <- read.table(bh_file, header=TRUE, sep="\t")
    bh <- bh[,1:2]
    colnames(bh) <- c("id", "active")
    
    fdr_inverse <- read.table(fdr_inverse_file, header=TRUE, sep="\t")
    ## collect value of qValue cutoff where samples were called present for qValue = 0.05 per individual sample
    FDRInverse_cutoff <- max(fdr_inverse$minimum_qValue[which(fdr_inverse$call == "present")])
    fdr_inverse <- fdr_inverse[,1:2]
    colnames(fdr_inverse) <- c("id", "active")
    
    zigzag <- read.table(zigzag_file, header=TRUE, sep="\t")
    zigzag <- zigzag[,1:2]
    colnames(zigzag) <- c("id","active")
    
    selectTableBH <- merge(TP_TN, bh, by ="id")
    selectTableBH <- selectTableBH[order(selectTableBH$active) , ]
    
    selectTableFDR <- merge(TP_TN, fdr_inverse, by ="id")
    selectTableFDR <- selectTableFDR[order(selectTableFDR$active) , ]
    
    selectTableZigzag <- merge(TP_TN, zigzag, by ="id")
    selectTableZigzag <- selectTableZigzag[order(selectTableZigzag$active) , ]
    
    return(list(FDRInverse_cutoff, selectTableBH, selectTableFDR, selectTableZigzag))
    
  }
}
