source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

## save output
outputFolder <- file.path("./figures/Suppl_figure_5/")
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

### results TPM >= 2
tpm_threshold <- dplyr::filter(tpm_threshold, libraryId %in% librariesUsed)
tpm_threshold <- tpm_threshold %>% dplyr::select(libraryId, proportionCodingPresent, speciesId, organism)
tpm_threshold$approach <- "TPM threshold"
tpm_threshold$cutoff <- "TPM >= 2"

## results with deconvolution and using different p-values cut-off
deconv_pValue_0.001 <- dplyr::filter(deconv_pValue, cutoff == 0.001 & libraryId %in% librariesUsed)
deconv_pValue_0.001 <- merge(deconv_pValue_0.001, sampleInfo, by="libraryId")
deconv_pValue_0.001 <- deconv_pValue_0.001 %>% dplyr::select(libraryId, proportionCodingPresent, speciesId, approach, cutoff, organism)

deconv_pValue_0.01 <- dplyr::filter(deconv_pValue, cutoff == 0.01 & libraryId %in% librariesUsed)
deconv_pValue_0.01 <- merge(deconv_pValue_0.01, sampleInfo, by="libraryId")
deconv_pValue_0.01 <- deconv_pValue_0.01 %>% dplyr::select(libraryId, proportionCodingPresent, speciesId, approach, cutoff, organism)

deconv_pValue_0.05 <- dplyr::filter(deconv_pValue, cutoff == 0.05 & libraryId %in% librariesUsed)
deconv_pValue_0.05 <- merge(deconv_pValue_0.05, sampleInfo, by="libraryId")
deconv_pValue_0.05 <- deconv_pValue_0.05 %>% dplyr::select(libraryId, proportionCodingPresent, speciesId, approach, cutoff, organism)

## results without decovolution
without_deconv <- dplyr::filter(without_deconv, libraryId %in% librariesUsed)
without_deconv <- merge(without_deconv, sampleInfo, by="libraryId")
without_deconv <- without_deconv %>% dplyr::select(libraryId, proportionCodingPresent,speciesId, organism)
without_deconv$approach <- "Without deconvolution"
without_deconv$cutoff <- "0.05"

################ re-organize big table
allInfo <- rbind(tpm_threshold, deconv_pValue_0.001, deconv_pValue_0.01, deconv_pValue_0.05, without_deconv)
allInfo$approach <- gsub("pValueCutoff", "Deconvolution p-value approach", allInfo$approach)
allInfo$cutoff <- gsub("0.001", "p-value <= 0.001", allInfo$cutoff)
allInfo$cutoff <- gsub("0.01", "p-value <= 0.01", allInfo$cutoff)
allInfo$cutoff <- gsub("0.05", "p-value <= 0.05", allInfo$cutoff)

### Plot for all species
for (i in unique(allInfo$speciesId)) {
  
  speciesData <- dplyr::filter(allInfo, speciesId == i)
  numberLib <- round(nrow(speciesData)/5, 0)  ### approaches
  speciesName <- as.character(unique(allInfo$organism[allInfo$speciesId == i]))
  
  graph <- ggplot(speciesData, aes(x=approach, y=proportionCodingPresent, fill=cutoff)) +
    geom_boxplot() + ylim(0,100) + ggtitle(paste0(speciesName))+
    xlab(" ") + ylab(" % Protein coding genes present")+
    annotate("text", x=c(1,2,3), y=100, label= paste0(numberLib))+ 
    theme(panel.background = element_rect(fill = "white", colour = "black"))
  
  pdf(file = file.path(outputFolder, paste0("Species_", speciesName, ".pdf")),width = 10, height = 4) 
  print(graph)
  dev.off()
  
}

## save output for stats
outputFolderStats <- file.path("./stats_info/All_species/")
if (!dir.exists(outputFolderStats)){
  dir.create(outputFolderStats)
} else {
  print("Already exists!")
}

## export table with information about % coding genes per library per method and cutoff
write.table(allInfo, file = file.path(outputFolderStats, "CodingGenes_callPresent_per_Library.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
