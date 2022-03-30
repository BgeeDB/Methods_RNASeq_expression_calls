source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- file.path("./figures/Figure_3/")
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

### results TPM >= 2
tpm_threshold <- dplyr::filter(tpm_threshold, libraryId %in% librariesUsed)  %>% dplyr::select(libraryId, proportionCodingPresent, speciesId, organism)
tpm_threshold$approach <- "TPM threshold"
tpm_threshold$cutoff <- "TPM >= 2"

## results with deconvolution and using different p-values cut-off
deconv_pValue_0.001 <- dplyr::filter(deconv_pValue, cutoff == 0.001 & libraryId %in% librariesUsed)
deconv_pValue_0.001 <- merge(deconv_pValue_0.001, sampleInfo, by="libraryId") %>% dplyr::select(libraryId, proportionCodingPresent, speciesId, approach, cutoff, organism)

deconv_pValue_0.01 <- dplyr::filter(deconv_pValue, cutoff == 0.01 & libraryId %in% librariesUsed)
deconv_pValue_0.01 <- merge(deconv_pValue_0.01, sampleInfo, by="libraryId")  %>% dplyr::select(libraryId, proportionCodingPresent, speciesId, approach, cutoff, organism)

deconv_pValue_0.05 <- dplyr::filter(deconv_pValue, cutoff == 0.05 & libraryId %in% librariesUsed)
deconv_pValue_0.05 <- merge(deconv_pValue_0.05, sampleInfo, by="libraryId")  %>% dplyr::select(libraryId, proportionCodingPresent, speciesId, approach, cutoff, organism)

## results without decovolution
without_deconv <- dplyr::filter(without_deconv, libraryId %in% librariesUsed)
without_deconv <- merge(without_deconv, sampleInfo, by="libraryId")  %>% dplyr::select(libraryId, proportionCodingPresent,speciesId, organism)
without_deconv$approach <- "Without deconvolution"
without_deconv$cutoff <- "0.05"

################ re-organize big table and plot just esox lucius species = 8010
allInfo <- rbind(tpm_threshold, deconv_pValue_0.001, deconv_pValue_0.01, deconv_pValue_0.05, without_deconv)
allInfo$approach <- gsub("pValueCutoff", "Deconvolution p-value approach", allInfo$approach)
allInfo$cutoff <- gsub("0.001", "p-value <= 0.001", allInfo$cutoff)
allInfo$cutoff <- gsub("0.01", "p-value <= 0.01", allInfo$cutoff)
allInfo$cutoff <- gsub("0.05", "p-value <= 0.05", allInfo$cutoff)

speciesData <- dplyr::filter(allInfo, speciesId == "8010")
numberLib <- round(nrow(speciesData)/5, 0)  ### approaches
speciesName <- as.character(unique(allInfo$organism[allInfo$speciesId == "8010"]))
  
graph <- ggplot(speciesData, aes(x=approach, y=proportionCodingPresent, fill=cutoff)) +
    geom_boxplot() + ylim(0,100) + ggtitle(paste0(speciesName))+
    xlab(" ") + ylab(" % Protein coding genes present")+
    annotate("text", x=c(1,2,3), y=100, label= paste0(numberLib))+ 
    theme(panel.background = element_rect(fill = "white", colour = "black"))
  
pdf(file = file.path(outputFolder, paste0("Protein_Coding_Present_", speciesName, ".pdf")),width = 10, height = 4) 
print(graph)
dev.off()
