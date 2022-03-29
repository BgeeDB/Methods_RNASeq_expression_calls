source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

## select just Mouse and Human data
species <- c(10090, 9606)
RNA_scRNA <- dplyr::filter(RNA_scRNA, speciesId %in% species)
RNA_scRNA$organism <- ifelse(RNA_scRNA$speciesId == "10090", "Mus musculus", "Homo sapiens")

outputFolder <- "./figures/Figure_6/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

pdf(file=file.path(outputFolder, "Compare_approach_rna_vs_scRNA_human_Mouse.pdf"),width=10, height=6)
g2 <- ggplot(RNA_scRNA, aes(organism, proportionCodingPresent, fill = dataType)) + 
  geom_boxplot() + ylim(0,100) + xlab(" ") + ylab("% Protein coding genes present") + ggtitle(" ")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))
g2
dev.off()
