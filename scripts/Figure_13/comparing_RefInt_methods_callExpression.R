source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- "./figures/Figure_13/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

callExpression <- merge(callExpression, rnaSeqLib, by ="libraryId")
comparingCalls_0.05 <- dplyr::filter(callExpression, cutoff == 0.05)
comparingCalls_0.05$organism <- NA
comparingCalls_0.05$organism <- ifelse(comparingCalls_0.05$speciesId == "7227", "Drosophila melanogaster", "NA")
comparingCalls_0.05$organism <- ifelse(comparingCalls_0.05$speciesId == "10090", "Mus musculus", comparingCalls_0.05$organism)
comparingCalls_0.05$organism <- ifelse(comparingCalls_0.05$speciesId == "9606", "Homo sapiens", comparingCalls_0.05$organism)

pdf(file = file.path(outputFolder,"Calls_pValue_vs_qValue.pdf"),width = 10, height = 4) 
g1 <- ggplot(comparingCalls_0.05, aes(x=as.factor(organism), y=proportionCodingPresent, fill=approach)) + 
  geom_boxplot() + ylim(0,100) + ggtitle("Calls with p-value versus q-value approach") +
  xlab(" ") + ylab("% Protein coding genes")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"))
g1
dev.off()
