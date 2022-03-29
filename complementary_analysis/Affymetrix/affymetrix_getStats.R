## Script used to provide information about call of expression in Mouse data for Affymetrix (table 1 of the paper)

outputFolder <- "./stats_info/Affymetrix_stats/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

musMusculus <- Bgee$new(species = "Mus_musculus", dataType = "affymetrix", release = "15_0")
musMusculusAnn <- getAnnotation(musMusculus)
## gene to biotype info
gene2biotype <- read.table(mouse_gene2biotype, header=TRUE, sep="\t")
colnames(gene2biotype)[1]  <- "Gene.ID"

## experiments
expMusmusculus <- read.table(paste0(here(),"/Mus_musculus_Bgee_15_0/Mus_musculus_Affymetrix_experiments.tsv"), header=TRUE, sep="\t")
dim(expMusmusculus)
libMusmusculus <- read.table(paste0(here(),"/Mus_musculus_Bgee_15_0/Mus_musculus_Affymetrix_chips.tsv"), header=TRUE, sep="\t")
numberSamples_gcRMA <- length(libMusmusculus$Chip.ID[libMusmusculus$Normalization.type == "gcRMA"])
numberSamples_mas5 <- length(libMusmusculus$Chip.ID[libMusmusculus$Normalization.type == "MAS5"])

collectChip_gcRMA <-libMusmusculus$Chip.ID[libMusmusculus$Normalization.type == "gcRMA"]
collectChip_MAS5 <- libMusmusculus$Chip.ID[libMusmusculus$Normalization.type == "MAS5"]

## getData for mouse
options(timeout=500)
musMusculusData <- getData(musMusculus)

gcRMAData <- dplyr::filter(musMusculusData, Chip.ID %in% collectChip_gcRMA)
mas5Data <- dplyr::filter(musMusculusData, Chip.ID %in% collectChip_MAS5)

## collect stats for MAS5 chip.ID
collectInfoMas5 <- c()
for (i in unique(mas5Data$Chip.ID)) {
  
  libraryID <- dplyr::filter(mas5Data, Chip.ID == i)
  libraryID <- merge(libraryID, gene2biotype, by="Gene.ID")
  libInfo <- (nrow(dplyr::filter(libraryID, biotype == "protein_coding", Detection.flag == "present" | Detection.flag == "marginal"))/nrow(dplyr::filter(libraryID, biotype == "protein_coding")))*100
  collectLib <- cbind(i,libInfo)
  collectInfoMas5 <- rbind(collectInfoMas5, collectLib)
}
colnames(collectInfoMas5) <- c("Chip.ID", "Proportion_Coding_present_marginal")
collectInfoMas5 <- as.data.frame(collectInfoMas5)
collectInfoMas5$Proportion_Coding_present_marginal <- as.numeric(collectInfoMas5$Proportion_Coding_present_marginal)

write.table(collectInfoMas5, file = file.path(outputFolder,"Affy_mas5_Proportion_ProteinCoding_per_chipID.tsv"), sep="\t", quote = FALSE, row.names = FALSE)
getInfo_mas5 <- summary(collectInfoMas5$Proportion_Coding_present_marginal)
capture.output(getInfo_mas5, file = file.path(outputFolder,"Affy_mas5_summary.tsv"))


## collect stats for gcRMA chip.ID
collectInfogcRMA <- c()
for (i in unique(gcRMAData$Chip.ID)) {
  
  libraryID <- dplyr::filter(gcRMAData, Chip.ID == i)
  libraryID <- merge(libraryID, gene2biotype, by="Gene.ID")
  libInfo <- (nrow(dplyr::filter(libraryID, biotype == "protein_coding", Detection.flag == "present" | Detection.flag == "marginal"))/nrow(dplyr::filter(libraryID, biotype == "protein_coding")))*100
  collectLib <- cbind(i,libInfo)
  collectInfogcRMA <- rbind(collectInfogcRMA, collectLib)
  
}
colnames(collectInfogcRMA) <- c("Chip.ID", "Proportion_Coding_present_marginal")
collectInfogcRMA <- as.data.frame(collectInfogcRMA)
collectInfogcRMA$Proportion_Coding_present_marginal <- as.numeric(collectInfogcRMA$Proportion_Coding_present_marginal)
write.table(collectInfogcRMA, file = file.path(outputFolder,"Affy_gcRMA_Proportion_ProteinCoding_per_chipID.tsv"), sep="\t", quote = FALSE, row.names = FALSE)
getInfo_gcRMA <- summary(collectInfogcRMA$Proportion_Coding_present_marginal)
capture.output(getInfo_gcRMA, file = file.path(outputFolder,"Affy_gcRMA_summary.tsv"))

