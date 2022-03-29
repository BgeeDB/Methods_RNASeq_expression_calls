source("./scripts//checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")
source("./scripts/InputFiles_for_benchmark.R")


## cutoff applied
cutoffRefInt_AllInt <- c(0.05, 0.01, 0.001, 0.00001)
cutoff_TPM <- c(0.5, 1, 2, 5, 10)

## libraries where analysis will be done
librarySamplesLiver <- c("ERX012346","ERX012361","ERX012364","ERX012370","ERX012375","ERX012377","ERX016336","ERX016337","ERX016340","ERX1174307","ERX1174311","ERX1174315","ERX1174323","ERX1174327",
                         "ERX2187489","ERX2187490","ERX2187491","SRX080232","SRX081917","SRX081918","SRX081919","SRX191151","SRX196268","SRX196276","SRX196285","SRX211608","SRX211609","SRX211610",
                         "SRX219281", "SRX262967","SRX262972","SRX1038930","SRX1603066","SRX1603067","SRX1603116","SRX1603117","SRX1603153","SRX1603154","SRX1603156","SRX1603157","SRX1603263",
                         "SRX1603264","SRX2751106", "SRX2751107","SRX2751108","SRX2751109","SRX2751110","SRX2751111","SRX2751112","SRX2751113","SRX2751114","SRX2751115","SRX2751116","SRX2751117",
                         "SRX2751130","SRX2751131","SRX2751132", "SRX2751133","SRX2751134","SRX2751135","SRX2751136","SRX2751137","SRX2751138","SRX2751139","SRX2751140","SRX2751141","SRX3050418",
                         "SRX3050419","SRX3050420","SRX3050421","SRX3050422","SRX3050423", "SRX3050424","SRX3050425","SRX3050426","SRX3050427","SRX3050428","SRX3050429","SRX3050430","SRX3050431",
                         "SRX3050432","SRX3050433","SRX3050434","SRX3050435","SRX3050436","SRX3050437","SRX3050438","SRX3050439","SRX3050440")

## save output
outputFolder <- file.path("./figures/Figure_2/")
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

#########################################################################################################3
## Function to calculate TPR and FPR per library
collectInfoMethods <- function(TP_TN_file, folder, libraryId, methodApplied, cutoff){
  
  TP_TN <- read.table(TP_TN_file, header=TRUE, sep="\t")
  colnames(TP_TN) <- c("id","callMarkers")
  
  collectInfo <- c()
  for (lib in libraryId) {
    print(lib)
    cutoff <- as.numeric(cutoff)
    
    if (methodApplied == "Ref_intergenic" | methodApplied == "TPM"){
      readFile <- read.table(file.path(folder, lib, "gene_level_abundance+calls_0.05.tsv"), header=TRUE, sep="\t")
    } else if (methodApplied == "All_intergenic") {
      readFile <- read.table(file.path(folder, lib, "gene_level_abundance+calls.tsv"), header=TRUE, sep="\t")
    } else {
      message("Method not recognized! You should specify: Ref_intergnic, TPM or All_intergenic.")
    }
    
    if (methodApplied == "Ref_intergenic" | methodApplied == "All_intergenic"){
      readFile$call <- ifelse(readFile$pValue <= cutoff, "present", "absent")
    } else {
      readFile$call <- ifelse(readFile$abundance >= cutoff, "present", "absent")
    }
    readFile$call <- readFile$call %>% replace_na('absent')
    selectTable <- merge(TP_TN, readFile, by ="id")
    
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
      collectInfo <- rbind(collectInfo,TP_FP_Info)
  }
    return(collectInfo)
}

## perform the analysis
allData_RefInt_allInt <- c()
for (i in cutoffRefInt_AllInt) {
  print(i)
  getanalysis_AllInt <- collectInfoMethods(TP_TN_file = TP_TN_file, folder = folderFiles_AllIntergenic, libraryId = librarySamplesLiver, methodApplied = "All_intergenic", cutoff = i)
  getanalysis_AllInt$Approach <- "Without_deconvolution"
  getanalysis_AllInt$cutoff <- i
  
  getanalysis_RefInt <- collectInfoMethods(TP_TN_file = TP_TN_file, folder = folderFiles_RefIntergenic, libraryId = librarySamplesLiver, methodApplied = "Ref_intergenic", cutoff = i)
  getanalysis_RefInt$Approach <- "Deconvolution"
  getanalysis_RefInt$cutoff <- i
  
  getAnalysis_forCutoff <- rbind(getanalysis_AllInt, getanalysis_RefInt)
  allData_RefInt_allInt <- rbind(allData_RefInt_allInt, getAnalysis_forCutoff)
}

allData_TPM <- c()
for (i in cutoff_TPM) {
  print(i)
  getanalysis_TPM <- collectInfoMethods(TP_TN_file = TP_TN_file, folder = folderFiles_RefIntergenic, libraryId = librarySamplesLiver, methodApplied = "TPM", cutoff = i)
  getanalysis_TPM$Approach <- "TPM"
  getanalysis_TPM$cutoff <- i
  
  allData_TPM <- rbind(allData_TPM, getanalysis_TPM)
}

## plot the TPR vs. FPR
pdf(file=file.path(outputFolder, "TRP_FPR_all_Mouse_liver_samples.pdf"),width=10, height=6)
allData_RefInt_allInt <- reshape2::melt(allData_RefInt_allInt, id.vars = c("Approach", "cutoff"),
                    measure.vars = c("TPR", "FPR"))

graphic_A <- ggplot(allData_RefInt_allInt, aes(x = Approach, y = value, linetype = factor(variable), fill = factor(cutoff))) + 
  geom_boxplot()+
  scale_linetype_manual(name = "Rates", values = c("solid", "dotted"))+
  guides(fill=guide_legend("cutoff")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))+ ggtitle("Intergenic approaches")

allData_TPM <- reshape2::melt(allData_TPM, id.vars = c("Approach", "cutoff"),
                              measure.vars = c("TPR", "FPR"))

graphic_B <- ggplot(allData_TPM, aes(x = Approach, y = value, linetype = factor(variable), fill = factor(cutoff))) + 
  geom_boxplot()+
  scale_linetype_manual(name = "Rates", values = c("solid", "dotted"))+
  guides(fill=guide_legend("cutoff")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))+ ggtitle("TPM approach")

grid.arrange(graphic_A, graphic_B, nrow = 1)
dev.off()
 




