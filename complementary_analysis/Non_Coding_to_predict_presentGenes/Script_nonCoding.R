source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- file.path("./stats_info/NonCoding_predict_genePresent/")
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

tpm2 <- dplyr::filter(tpm_threshold, libraryId %in% all_lib_bgee15$libraryId)

intDev <- dplyr::filter(deconv_pValue, cutoff == "0.05")
intDev <- dplyr::filter(intDev, libraryId %in% all_lib_bgee15$libraryId)

intNODev <- dplyr::filter(without_deconv, libraryId %in% all_lib_bgee15$libraryId)

uniqueSpecies <- unique(tpm2$organism)

pathData <- "./data/Species_information/"
listFiles <- list.files(pathData, pattern = "*.gene2biotype", full.names = TRUE)

getProportionNoCodingGenes <- c()
for (file in listFiles) {
  
  nameFile <- basename(file)
  print(nameFile)
  nameSpecies <- str_extract(nameFile, "[^.]+")
  nameSpecies <- gsub("_", " ", nameSpecies)
  
  readFile <- read.table(file, header = TRUE, sep="\t")
  readFile <- na.omit(readFile) 
  non_coding_annot <- nrow(dplyr::filter(readFile, biotype != "protein_coding"))
  proportioNonCodingGenome <- (non_coding_annot / nrow(readFile))*100
  print(proportioNonCodingGenome)
  
  infoSpecies <- c(nameSpecies, non_coding_annot, proportioNonCodingGenome)
  getProportionNoCodingGenes <- rbind(getProportionNoCodingGenes, infoSpecies)
}
getProportionNoCodingGenes <- as.data.frame(getProportionNoCodingGenes)  
colnames(getProportionNoCodingGenes) <- c("organism", "Number_NonCoding", "ProportionNonCoding")

tpm2 <- merge(tpm2, getProportionNoCodingGenes, by="organism")
tpm2$ProportionNonCoding <-as.numeric(tpm2$ProportionNonCoding)

getLibSpeciesInfo <- data.frame(tpm2$libraryId, tpm2$organism, tpm2$Number_NonCoding, tpm2$ProportionNonCoding)
colnames(getLibSpeciesInfo) <- c("libraryId", "organism", "Number_NonCoding", "ProportionNonCoding")

intDev <- merge(intDev, getLibSpeciesInfo, by = "libraryId")
intDev$ProportionNonCoding <- as.numeric(intDev$ProportionNonCoding)

intNODev <- merge(intNODev, getLibSpeciesInfo, by = "libraryId")
intNODev$ProportionNonCoding <- as.numeric(intNODev$ProportionNonCoding)


### get median of protein coding genes present per species
collectInfoMedian <- c()
for (species in uniqueSpecies) {
  
  allCodingTPM <- dplyr::filter(tpm2, organism == species)
  medianCodingTPM <- median(allCodingTPM$proportionCodingPresent)
  
  allCodingintDev <- dplyr::filter(intDev, organism == species)
  medianCodingIntDev <- median(allCodingintDev$proportionCodingPresent)
  
  allCodingWithoutintDev <- dplyr::filter(intNODev, organism == species)
  medianCodingWithoutintDev <- median(allCodingWithoutintDev$proportionCodingPresent)
  
  allMethods <- c(species, medianCodingTPM, medianCodingIntDev, medianCodingWithoutintDev, unique(allCodingintDev$ProportionNonCoding))
  
  collectInfoMedian <- rbind(collectInfoMedian, allMethods)
}
colnames(collectInfoMedian) <- c("Species", "medianCodingTPM", "medianCodingIntDev", "medianCodingWithoutintDev", "Proportion_NoCoding")
collectInfoMedian <- as.data.frame(collectInfoMedian)
collectInfoMedian$medianCodingTPM <- as.numeric(collectInfoMedian$medianCodingTPM)
collectInfoMedian$medianCodingIntDev <- as.numeric(collectInfoMedian$medianCodingIntDev)
collectInfoMedian$medianCodingWithoutintDev <- as.numeric(collectInfoMedian$medianCodingWithoutintDev)
collectInfoMedian$Proportion_NoCoding <- as.numeric(collectInfoMedian$Proportion_NoCoding)

## correlation using all samples/libraries
tpm_cor_pearson <- cor.test(collectInfoMedian$medianCodingTPM, collectInfoMedian$Proportion_NoCoding, method = "pearson")
tpm_cor_spearman <- cor.test(collectInfoMedian$medianCodingTPM, collectInfoMedian$Proportion_NoCoding, method = "spearman")

intDec_cor_pearson <- cor.test(collectInfoMedian$medianCodingIntDev, collectInfoMedian$Proportion_NoCoding, method = "pearson")
intDec_cor_spearman <- cor.test(collectInfoMedian$medianCodingIntDev, collectInfoMedian$Proportion_NoCoding, method = "spearman")

withoutDec_cor_pearson <- cor.test(collectInfoMedian$medianCodingWithoutintDev, collectInfoMedian$Proportion_NoCoding, method = "pearson")
withoutDec_cor_spearman <- cor.test(collectInfoMedian$medianCodingWithoutintDev, collectInfoMedian$Proportion_NoCoding, method = "spearman")

collectAllStats <- list(tpm_cor_pearson, tpm_cor_spearman, intDec_cor_pearson, intDec_cor_spearman, withoutDec_cor_pearson, withoutDec_cor_spearman)
names(collectAllStats) <- c("Pearson correlation: median coding TPM vs. Proportion non coding", "Spearman correlation: median coding TPM vs. Proportion non coding", 
                            "Pearson correlation: median coding deconvolution vs. Proportion non coding", "Spearman correlation: median coding deconvolution vs. Proportion non coding",
                            "Pearson correlation: median coding without deconvolution vs. Proportion non coding", "Spearman correlation: median coding without deconvolution vs. Proportion non coding")
capture.output(collectAllStats, file = file.path(outputFolder,"NonCoding_to_predictPresentGenes.tsv"))
