source("./scripts/checkLibraries.R")

outputFolder <- "./stats_info/Blood_stats/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

### blood samples from GTEX
bloodSamples <- fread("./complementary_analysis/4-blood_human_data/Blood_samples_Bgee15.tsv", header=TRUE, sep="\t")

## select just GTEx samples
gtexSamples <- dplyr::filter(bloodSamples, Experiment.ID == "SRP012682")

collectInformation <- c()
for (i in unique(gtexSamples$Library.ID)) {
  
  genesBlood <- c("ENSG00000206172", "ENSG00000188536", "ENSG00000244734")
  
  sumReads <- sum(gtexSamples$Read.count[gtexSamples$Library.ID == i])
  genesExpressed <- nrow(dplyr::filter(gtexSamples, Library.ID == i & Detection.flag == "present"))
  
  ## select target blood genes
  bloodGenes <- dplyr::filter(gtexSamples, Library.ID == i & Gene.ID %in% genesBlood)
  sumReadsBlood <- sum(bloodGenes$Read.count)
  
  getProportion <- (sumReadsBlood/sumReads)*100
  
  collectLibraryInfo <- cbind(i, getProportion, genesExpressed, sumReadsBlood)
  collectInformation <- rbind(collectInformation, collectLibraryInfo)
}
collectInformation <- as.data.frame(collectInformation)
collectInformation$getProportion <- as.numeric(collectInformation$getProportion)
collectInformation$genesExpressed <- as.numeric(collectInformation$genesExpressed)
collectInformation$sumReadsBlood <- as.numeric(collectInformation$sumReadsBlood)
colnames(collectInformation) <- c("libraryId", "getProportion_bloodGenes", "Genes_expressed_library", "reads_blood_genes")


summaryBlood <- summary(collectInformation$getProportion_bloodGenes)
## correlation between HBs proportion and genes called present
getCorrelation <- cor(collectInformation$reads_blood_genes, collectInformation$Genes_expressed_library)

collectAllStats <- list(summaryBlood, getCorrelation)
names(collectAllStats) <- c("Summary stats for proportion of blood genes", "Correlation Sum reads for blood genes versus expressed genes")
capture.output(collectAllStats, file = file.path(outputFolder,"Stats_blood_samples_from_GTEx.tsv"))
