## This script is used to plot the proportion of non coding genes versus the proportion of selected reference intergenic regions per species.
source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- "./figures/Suppl_figure_11/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

pathData <- "./data/Species_information/"
listFiles <- list.files(pathData, pattern = "*.gene2biotype", full.names = TRUE)

#### GET NON CODING 
getProportionNoCodingGenes <- c()
for (file in listFiles) {
  
  nameFile <- basename(file)
  nameSpecies <- str_extract(nameFile, "[^.]+")
  nameSpecies <- gsub("_", " ", nameSpecies)
  
  readFile <- read.table(file, header = TRUE, sep="\t")
  readFile <- na.omit(readFile) 
  non_coding_annot <- nrow(dplyr::filter(readFile, biotype != "protein_coding"))
  proportioNonCodingGenome <- (non_coding_annot / nrow(readFile))*100
  
  infoSpecies <- c(nameSpecies, non_coding_annot, proportioNonCodingGenome)
  getProportionNoCodingGenes <- rbind(getProportionNoCodingGenes, infoSpecies)
}
getProportionNoCodingGenes <- as.data.frame(getProportionNoCodingGenes)  
colnames(getProportionNoCodingGenes) <- c("Species_name", "Number_NonCoding", "Proportion_NonCoding")

getProportionNoCodingGenes <- getProportionNoCodingGenes[order(getProportionNoCodingGenes$Species_name),]
speciesID <- speciesInfoCollected[order(speciesInfoCollected$Species_name),]
getProportionNoCodingGenes$Species_id <- speciesID$Species_id

## file per species with Ref.Int + Other intergenic
intergenic <- "./data/RefInt+OtherInt_All_species/"
species <- c(9606, 10090, 7955, 7227, 6239, 9598, 9597, 9593, 9544, 10116, 9913, 9823, 9796, 9986, 9615, 9685, 
             10141, 13616, 9258, 9031, 28377, 8364, 7237, 7240, 7740,7897, 7918, 7936, 7994, 8010, 8030, 8049, 
             8081, 8090, 8154, 8355, 9103, 9483, 9531, 9541, 9545, 9555, 9925, 9940, 9974, 10181, 30608, 32507, 52904, 60711, 69293, 105023)

## for all species collect information
allSpeciesInfo <- c()
## loop throught species
collectProportion <- c()
for (i in species) {
  
  readFile <- fread(paste0(intergenic, "All_regions_", i, ".tsv"), header=TRUE, sep="\t")
  
  ## Number of other regions
  allOtherIntergenic <- nrow(dplyr::filter(readFile, markName == "Other Ref. Intergenic"))
  ## Number of ref intergenic regions
  allRefIntergenic <- nrow(dplyr::filter(readFile, markName == "Ref. Intergenic"))
  
  ## all intergenic regions (other + ref)
  allIntergenic <- allOtherIntergenic+allRefIntergenic
  proportionInSpecies_otherInt <- (allOtherIntergenic/allIntergenic)*100
  proportionInSpecies_refInt <- (allRefIntergenic/allIntergenic)*100
  
  collectInfoSpecies <- c(i, proportionInSpecies_otherInt, proportionInSpecies_refInt)
  
  collectInfoChro <- c()
  for (j in unique(readFile$chr)) {
    
    ## collect info per chr
    refInt_chr <- nrow(dplyr::filter(readFile, chr == j & markName == "Ref. Intergenic"))
    proportionRefInt <- (refInt_chr/allRefIntergenic)*100
    RefBind <- cbind(refInt_chr, proportionRefInt)
    
    otherInt_chr <- nrow(dplyr::filter(readFile, chr == j & markName == "Other Ref. Intergenic"))
    proportionOtherInt <- (otherInt_chr/allOtherIntergenic)*100
    OtherBind <- cbind(otherInt_chr, proportionOtherInt)
    
    dataInfo <- rbind(RefBind, OtherBind)
    row.names(dataInfo) <- c("RefBind", "OtherBind")
    colnames(dataInfo) <- c("Number_regions", "proportion_regions")
    dataInfo <- as.data.frame(dataInfo)
    dataInfo$chr <- paste0(j)
    dataInfo$Region <- c("Ref_intergenic", "Other_intergenic")
    dataInfo$speciesId <- i
    
    collectInfoChro <- rbind(collectInfoChro, dataInfo)
  }
  allSpeciesInfo <- rbind(allSpeciesInfo, collectInfoChro)
  collectProportion <- rbind(collectProportion, collectInfoSpecies)
}
colnames(collectProportion) <- c("Species_id", "Proportion_OtherIntergenic", "Proportion_RefIntergenic")
collectProportion <- as.data.frame(collectProportion)

allInfo <- merge(getProportionNoCodingGenes, collectProportion, by="Species_id")
allInfo$Proportion_NonCoding <- as.numeric(allInfo$Proportion_NonCoding)


pdf(file = file.path(outputFolder, "/ProportionNoncoding_vs_ProportionRefIntergenic.pdf"), width = 16, height = 7)
p <- ggplot(data=allInfo, aes(x=Proportion_NonCoding, y=Proportion_RefIntergenic, label = Species_id)) +
  geom_point(stat="identity")+ ylim(0,100)+ xlim(0,100) +
  geom_text(hjust=1.2, vjust=0)+
  xlab("Proportion of Non Coding genes") + ylab("Proportion of Reference Intergenic region selected") +
  ggtitle("Reference intergenic regions vs. Non coding genes") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),legend.position = "none")
p
dev.off()





