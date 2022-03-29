source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_for_benchmark.R")


librarySamplesDrosophilaTestis <- c("SRX1720958", "SRX1720957", "SRX493999", "SRX493950", "SRX109279", "SRX109278")

librarySamplesHumanLung <- c("SRX4794949","SRX4794942","SRX4794906","SRX4794874","SRX4794873","SRX4794872","SRX4794871",
                             "SRX4794870","SRX1020489","SRX1020488","SRX815518","SRX815517","SRX815516","SRX815515",
                             "SRX815514","SRX815513","ERX288640","ERX288620","ERX288615","ERX288597","ERX288596",
                             "ERX288531","ERX288496","ERX288491","ERX011227", "ERX011222")

## save output
outputFolder <- file.path("./figures/Suppl_figure_4/")
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

##################################################################################################################################################################################################################
### Collect input data and calculate TPR and FPR
collectInfoMethods <- function(TP_TN_file, folder, libraryId, methodApplied, cutoff){
  
  TP_TN <- read.table(TP_TN_file, header=TRUE, sep="\t")
  colnames(TP_TN) <- c("id","callMarkers")
  
  collectInfo <- c()
  for (lib in libraryId) {
    print(lib)
    cutoff <- as.numeric(cutoff)
    
    if (methodApplied == "Ref_intergenic"){
      readFile <- read.table(file.path(folder, lib, "gene_level_abundance+calls_0.05.tsv"), header=TRUE, sep="\t")
      readFile$calculateBH <- p.adjust(readFile$pValue, method = "BH") 
      readFile$call <- ifelse(readFile$pValue <= cutoff, "present", "absent")
      readFile$callBH <- ifelse(readFile$calculateBH <= cutoff, "present", "absent")
      readFile$call <- readFile$call %>% replace_na('absent')
      readFile$callBH <- readFile$callBH %>% replace_na('absent')
      selectTable <- merge(TP_TN, readFile, by ="id")
    } else if (methodApplied == "All_intergenic"){
      readFile <- read.table(file.path(folder, lib, "gene_level_abundance+calls.tsv"), header=TRUE, sep="\t")
      readFile$calculateBH <- p.adjust(readFile$pValue, method = "BH") 
      readFile$call <- ifelse(readFile$pValue <= cutoff, "present", "absent")
      readFile$callBH <- ifelse(readFile$calculateBH <= cutoff, "present", "absent")
      readFile$call <- readFile$call %>% replace_na('absent')
      readFile$callBH <- readFile$callBH %>% replace_na('absent')
      selectTable <- merge(TP_TN, readFile, by ="id")
    } else if (methodApplied == "TPM"){
      readFile <- read.table(file.path(folder, lib, "gene_level_abundance+calls_0.05.tsv"), header=TRUE, sep="\t")
      readFile$call <- ifelse(readFile$abundance >= cutoff, "present", "absent")
      readFile$call <- readFile$call %>% replace_na('absent')
      selectTable <- merge(TP_TN, readFile, by ="id")
    } else {
      message("Method applied not recognized!")
    }
    
    if (methodApplied == "Ref_intergenic" | methodApplied == "All_intergenic"){
      ## present in control markers file and also in the method
      TP_pvalue <- dplyr::filter(selectTable, callMarkers == "present" & call == "present")
      TP_pvalue <- nrow(TP_pvalue)
      TP_BH <- dplyr::filter(selectTable, callMarkers == "present" & callBH == "present")
      TP_BH <- nrow(TP_BH)
      ## absent in control markers file and present in method
      FP_pvalue <- dplyr::filter(selectTable, callMarkers == "absent" & call == "present")
      FP_pvalue <- nrow(FP_pvalue)
      FP_BH <- dplyr::filter(selectTable, callMarkers == "absent" & callBH == "present")
      FP_BH <- nrow(FP_BH)
      ## absent in the control markers file and absent in the method
      TN_pvalue <- dplyr::filter(selectTable, callMarkers == "absent" & call == "absent")
      TN_pvalue <- nrow(TN_pvalue)
      TN_BH <- dplyr::filter(selectTable, callMarkers == "absent" & callBH == "absent")
      TN_BH <- nrow(TN_BH)
      ## present in the control markers file and absent in the method
      FN_pvalue <- dplyr::filter(selectTable, callMarkers == "present" & call == "absent")
      FN_pvalue <- nrow(FN_pvalue)
      FN_BH <- dplyr::filter(selectTable, callMarkers == "present" & callBH == "absent")
      FN_BH <- nrow(FN_BH)
      ## Calculate ratios
      TPR_pValue <- TN_pvalue/(TN_pvalue + FN_pvalue)
      FPR_pvalue <- FP_pvalue/(FP_pvalue+TN_pvalue)
      TPR_BH <- TN_BH/( TN_BH + FN_BH)
      FPR_BH <- FP_BH/(FP_BH+TN_BH)
      
      TP_FP_Info <- data.frame(TP_pvalue,TP_BH, FP_pvalue, FP_BH, TN_pvalue, TN_BH, FN_pvalue, FN_BH, TPR_pValue, FPR_pvalue, TPR_BH, FPR_BH)
    } else {
      ## present in control markers file and also in the method
      TP <- dplyr::filter(selectTable, callMarkers == "present" & call == "present")
      TP <- nrow(TP)
      ## absent in control markers file and present in method
      FP <- dplyr::filter(selectTable, callMarkers == "absent" & call == "present")
      FP <- nrow(FP)
      ## absent in the control markers file and absent in the method
      TN <- dplyr::filter(selectTable, callMarkers == "absent" & call == "absent")
      TN <- nrow(TN)
      ## present in the control markers file and absent in the method
      FN <- dplyr::filter(selectTable, callMarkers == "present" & call == "absent")
      FN <- nrow(FN)
      ## Calculate ratios
      TPR <- TN/(TN + FN)
      FPR <- FP/(FP+TN)
      
      TP_FP_Info <- data.frame(TP, FP, TN, FN, TPR, FPR)
    }
    collectInfo <- rbind(collectInfo,TP_FP_Info)
  }
  collectInfo$approach <- ifelse(methodApplied == "All_intergenic", "Without_deconvolution" , ifelse(methodApplied == "Ref_intergenic", "Deconvolution", "TPM"))
  collectInfo$cutoff <- cutoff
  return(collectInfo)
}


################################# Plot supplementary figures: with and without BH correction for Human lung #############################
##### Collect stats to without deconvolution for different pValues cut-off
cutoff <- c(0.05, 0.01, 0.001, 0.00001)
collectAllInt <- c()
for (i in cutoff) {
  getInfo <- collectInfoMethods(TP_TN_file = TP_TN_human_lung, folder = folderHumanLungAllInt, libraryId = librarySamplesHumanLung, methodApplied = "All_intergenic", cutoff = i)
  collectAllInt <- rbind(collectAllInt, getInfo)
}

##### Collect stats to deconvolution method for different pValues cut-off
cutoff <- c(0.05, 0.01, 0.001, 0.00001)
collectRefInt <- c()
for (i in cutoff) {
  getInfo <- collectInfoMethods(TP_TN_file = TP_TN_human_lung, folder = folderHumanLungRefInt, libraryId = librarySamplesHumanLung, methodApplied = "Ref_intergenic", cutoff = i)
  collectRefInt <- rbind(collectRefInt, getInfo)
}

##### Collect stats to TPM method for different TPM threshold cut-off
cutoff <- c(0.5, 1, 2, 5, 10)
collectTPM <- c()
for (i in cutoff) {
  getInfo <- collectInfoMethods(TP_TN_file = TP_TN_human_lung, folder = folderHumanLungRefInt, libraryId = librarySamplesHumanLung, methodApplied = "TPM", cutoff = i)
  collectTPM <- rbind(collectTPM, getInfo)
}


allInterg_RefInt <- rbind(collectAllInt, collectRefInt)
allInterg_RefInt <- allInterg_RefInt %>% select(TPR_pValue, FPR_pvalue, TPR_BH, FPR_BH, approach, cutoff)
colnames(allInterg_RefInt) <- c("TPR", "FPR", "TPR_BH", "FPR_BH", "Approach", "cutoff")
allInterg_RefInt = melt(allInterg_RefInt, id.vars = c("Approach", "cutoff"),
                        measure.vars = c("TPR", "FPR","TPR_BH", "FPR_BH"))

tpm <- collectTPM %>% select(TPR, FPR,approach, cutoff)
colnames(tpm) <- c("TPR", "FPR", "Approach", "cutoff")
tpm_approach = melt(tpm, id.vars = c("Approach", "cutoff"),
                    measure.vars = c("TPR", "FPR"))

g1 <- ggplot(allInterg_RefInt, aes(x = Approach, y = value, linetype = factor(variable), fill = factor(cutoff))) + 
  geom_boxplot()+
  scale_linetype_manual(name = "Rates", values = c("solid", "dotted", "dashed", "longdash"))+
  guides(fill=guide_legend("cutoff")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))+ ggtitle("Intergenic approaches")

g2 <- ggplot(tpm_approach, aes(x = Approach, y = value, linetype = factor(variable), fill = factor(cutoff))) + 
  geom_boxplot()+
  scale_linetype_manual(name = "Rates", values = c("solid", "dotted"))+
  guides(fill=guide_legend("cutoff")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))+ ggtitle("TPM approach")

pdf(file = file.path(outputFolder,"HumanLung_with_and_without_BH_correction.pdf"),width = 16, height = 8) 
grid.arrange(g1, g2, nrow = 1)
dev.off()



################################# Plot supplementary figures: with and without BH correction for Drosophila Testis #############################
##### Collect stats to without deconvolution for different pValues cut-off
cutoff <- c(0.05, 0.01, 0.001, 0.00001)
collectAllInt <- c()
for (i in cutoff) {
  getInfo <- collectInfoMethods(TP_TN_file = TP_TN_drosophila_testis, folder = folderDrosTestisAllInt, libraryId = librarySamplesDrosophilaTestis, methodApplied = "All_intergenic", cutoff = i)
  collectAllInt <- rbind(collectAllInt, getInfo)
}

##### Collect stats to deconvolution method for different pValues cut-off
cutoff <- c(0.05, 0.01, 0.001, 0.00001)
collectRefInt <- c()
for (i in cutoff) {
  getInfo <- collectInfoMethods(TP_TN_file = TP_TN_drosophila_testis, folder = folderDrosTestisRefInt, libraryId = librarySamplesDrosophilaTestis, methodApplied = "Ref_intergenic", cutoff = i)
  collectRefInt <- rbind(collectRefInt, getInfo)
}

##### Collect stats to TPM method for different TPM threshold cut-off
cutoff <- c(0.5, 1, 2, 5, 10)
collectTPM <- c()
for (i in cutoff) {
  getInfo <- collectInfoMethods(TP_TN_file = TP_TN_drosophila_testis, folder = folderDrosTestisRefInt, libraryId = librarySamplesDrosophilaTestis, methodApplied = "TPM", cutoff = i)
  collectTPM <- rbind(collectTPM, getInfo)
}


allInterg_RefInt <- rbind(collectAllInt, collectRefInt)
allInterg_RefInt <- allInterg_RefInt %>% select(TPR_pValue, FPR_pvalue, TPR_BH, FPR_BH, approach, cutoff)
colnames(allInterg_RefInt) <- c("TPR", "FPR", "TPR_BH", "FPR_BH", "Approach", "cutoff")
allInterg_RefInt = melt(allInterg_RefInt, id.vars = c("Approach", "cutoff"),
                        measure.vars = c("TPR", "FPR","TPR_BH", "FPR_BH"))

tpm <- collectTPM %>% select(TPR, FPR,approach, cutoff)
colnames(tpm) <- c("TPR", "FPR", "Approach", "cutoff")
tpm_approach = melt(tpm, id.vars = c("Approach", "cutoff"),
                    measure.vars = c("TPR", "FPR"))

g1 <- ggplot(allInterg_RefInt, aes(x = Approach, y = value, linetype = factor(variable), fill = factor(cutoff))) + 
  geom_boxplot()+
  scale_linetype_manual(name = "Rates", values = c("solid", "dotted", "dashed", "longdash"))+
  guides(fill=guide_legend("cutoff")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))+ ggtitle("Intergenic approaches")

g2 <- ggplot(tpm_approach, aes(x = Approach, y = value, linetype = factor(variable), fill = factor(cutoff))) + 
  geom_boxplot()+
  scale_linetype_manual(name = "Rates", values = c("solid", "dotted"))+
  guides(fill=guide_legend("cutoff")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))+ ggtitle("TPM approach")

pdf(file = file.path(outputFolder,"DrosophilaTestis_with_and_without_BH_correction.pdf"),width = 16, height = 8) 
grid.arrange(g1, g2, nrow = 1)
dev.off()


