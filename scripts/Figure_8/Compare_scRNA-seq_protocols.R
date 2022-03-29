source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- "./figures/Figure_8/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

species <- c(10090, 9606)
scRNA_FL <- dplyr::filter(RNA_scRNA, speciesId %in% species & dataType == "scRNA-Seq")
scRNA_FL <- scRNA_FL %>% dplyr::select(libraryId, proportionGenicPresent, proportionCodingPresent, proportionIntergenicPresent, meanIntergenic, sdIntergenic, cutoffTPM, speciesId)
colnames(scRNA_FL)[7] <- "cutoff"
scRNA_FL$protocol <- "full-length"

### get results from target based protocols
targetBased <- targetBased %>% dplyr::select(libraryId, Proportion_genic_present, Proportion_coding_present, Proportion_intergenic_present, meanRefIntergenic, sdRefIntergenic, CPM_Threshold, species)
colnames(targetBased) <- c("libraryId","proportionGenicPresent", "proportionCodingPresent", "proportionIntergenicPresent", "meanIntergenic", "sdIntergenic", "cutoff", "speciesId")
targetBased$protocol <- "target-based"

allInfo_singleCell <- rbind(scRNA_FL, targetBased)
allInfo_singleCell$organism <- ifelse(allInfo_singleCell$speciesId == 10090, "Mus musculus", "Homo sapiens")

pdf(file=file.path(outputFolder,"Compare_singleCell_protocols_calls_human_Mouse.pdf"),width=10, height=6)
scRNASeq_protocols <- ggplot(allInfo_singleCell, aes(organism, proportionCodingPresent, fill = protocol)) + 
  geom_boxplot() + ylim(0,100) + xlab(" ") + ylab("% Protein coding genes present") + ggtitle(" ")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))
scRNASeq_protocols
dev.off()
