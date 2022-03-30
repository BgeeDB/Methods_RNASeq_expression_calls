source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_for_benchmark.R")

librarySamplesMouseLiver <- c("ERX012346","ERX012361","ERX012364","ERX012370","ERX012375","ERX012377","ERX016336","ERX016337","ERX016340","ERX1174307","ERX1174311","ERX1174315","ERX1174323","ERX1174327",
                              "ERX2187489","ERX2187490","ERX2187491","SRX080232","SRX081917","SRX081918","SRX081919","SRX191151","SRX196268","SRX196276","SRX196285","SRX211608","SRX211609","SRX211610",
                              "SRX219281", "SRX262967","SRX262972","SRX1038930","SRX1603066","SRX1603067","SRX1603116","SRX1603117","SRX1603153","SRX1603154","SRX1603156","SRX1603157","SRX1603263",
                              "SRX1603264","SRX2751106", "SRX2751107","SRX2751108","SRX2751109","SRX2751110","SRX2751111","SRX2751112","SRX2751113","SRX2751114","SRX2751115","SRX2751116","SRX2751117",
                              "SRX2751130","SRX2751131","SRX2751132", "SRX2751133","SRX2751134","SRX2751135","SRX2751136","SRX2751137","SRX2751138","SRX2751139","SRX2751140","SRX2751141","SRX3050418",
                              "SRX3050419","SRX3050420","SRX3050421","SRX3050422","SRX3050423", "SRX3050424","SRX3050425","SRX3050426","SRX3050427","SRX3050428","SRX3050429","SRX3050430","SRX3050431",
                              "SRX3050432","SRX3050433","SRX3050434","SRX3050435","SRX3050436","SRX3050437","SRX3050438","SRX3050439","SRX3050440")

## save output
outputFolder <- file.path("./figures/Suppl_figure_3/")
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
      readFile$callBH <- ifelse(readFile$calculateBH <= cutoff, "present", "absent")
      readFile$callBH <- readFile$callBH %>% replace_na('absent')
      selectTable <- merge(TP_TN, readFile, by ="id")
    } else if (methodApplied == "All_intergenic"){
      readFile <- read.table(file.path(folder, lib, "gene_level_abundance+calls.tsv"), header=TRUE, sep="\t")
      readFile$calculateBH <- p.adjust(readFile$pValue, method = "BH") 
      readFile$callBH <- ifelse(readFile$calculateBH <= cutoff, "present", "absent")
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
      TP <- dplyr::filter(selectTable, callMarkers == "present" & callBH == "present")
      TP <- nrow(TP)
      ## absent in control markers file and present in method
      FP <- dplyr::filter(selectTable, callMarkers == "absent" & callBH == "present")
      FP <- nrow(FP)
      ## absent in the control markers file and absent in the method
      TN <- dplyr::filter(selectTable, callMarkers == "absent" & callBH == "absent")
      TN <- nrow(TN)
      ## present in the control markers file and absent in the method
      FN <- dplyr::filter(selectTable, callMarkers == "present" & callBH == "absent")
      FN <- nrow(FN)
      ## Calculate ratios
      TPR <- TN/( TN + FN)
      FPR <- FP/(FP+TN)
      
      TP_FP_Info <- data.frame(TP, FP, TN, FN, TPR, FPR)
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

### Collect all the stats uwing the different methods
# Collect stats to without deconvolution for different pValues cut-off
cutoff <- c(0.05, 0.01, 0.001, 0.00001)
collectAllInt <- c()
for (i in cutoff) {
  getInfo <- collectInfoMethods(TP_TN_file = TP_TN_file, folder = folderFiles_AllIntergenic, libraryId = librarySamplesMouseLiver, methodApplied = "All_intergenic", cutoff = i)
  collectAllInt <- rbind(collectAllInt, getInfo)
}

# Collect stats to deconvolution method for different pValues cut-off
cutoff <- c(0.05, 0.01, 0.001, 0.00001)
collectRefInt <- c()
for (i in cutoff) {
  getInfo <- collectInfoMethods(TP_TN_file = TP_TN_file, folder = folderFiles_RefIntergenic, libraryId = librarySamplesMouseLiver, methodApplied = "Ref_intergenic", cutoff = i)
  collectRefInt <- rbind(collectRefInt, getInfo)
}

# Collect stats to TPM method for different TPM threshold cut-off
cutoff <- c(0.5, 1, 2, 5, 10)
collectTPM <- c()
for (i in cutoff) {
  getInfo <- collectInfoMethods(TP_TN_file = TP_TN_file, folder = folderFiles_RefIntergenic, libraryId = librarySamplesMouseLiver, methodApplied = "TPM", cutoff = i)
  collectTPM <- rbind(collectTPM, getInfo)
}


allInterg_RefInt <- rbind(collectAllInt, collectRefInt)
allInterg_RefInt <- allInterg_RefInt %>% select(TPR, FPR,approach, cutoff)
colnames(allInterg_RefInt) <- c("TPR", "FPR", "Approach", "cutoff")
allInterg_RefInt = melt(allInterg_RefInt, id.vars = c("Approach", "cutoff"),
                        measure.vars = c("TPR", "FPR"))

tpm <- collectTPM %>% select(TPR, FPR,approach, cutoff)
colnames(tpm) <- c("TPR", "FPR", "Approach", "cutoff")
tpm_approach = melt(tpm, id.vars = c("Approach", "cutoff"),
                    measure.vars = c("TPR", "FPR"))

g1 <- ggplot(allInterg_RefInt, aes(x = Approach, y = value, linetype = factor(variable), fill = factor(cutoff))) + 
  geom_boxplot()+
  scale_linetype_manual(name = "Rates", values = c("solid", "dotted"))+
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

pdf(file = file.path(outputFolder,"MouseLiver_with_BH_correction.pdf"),width = 16, height = 8) 
grid.arrange(g1, g2, nrow = 1)
dev.off()

