source("./scripts/checkLibraries.R")

outputFolder <- "./figures/Figure_9/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

pantherResults <- read.table("./data/PANTHER_EnrichmentAnalysis/Panther_results_BgeeDecv_CPM.tsv", header=TRUE, sep="\t")

## Overrepresentation pathways using CPM approach
randomCutoff <- pantherResults[pantherResults$method == "CPM",]
cpm <- ggplot(randomCutoff, aes(x=reorder(PANTHER_Pathways, log(upload_1.fold.Enrichment.)), y=log(upload_1.fold.Enrichment.))) + 
  geom_bar(stat = "identity") +
  coord_flip() + ggtitle(expression("Overrepresentation pathways: CPM">=1)) +
  xlab(" ") + ylab("Log(Fold Enrichment)") 

## Overrepresentation pathways using deconvoluiton approach
deconvolutionCalls <- pantherResults[pantherResults$method == "Deconvolution",]
deconvolution <- ggplot(deconvolutionCalls, aes(x=reorder(PANTHER_Pathways, log(upload_1.fold.Enrichment.)), y=log(upload_1.fold.Enrichment.))) + 
  geom_bar(stat = "identity") +
  coord_flip() + ggtitle("Overrepresentation pathways: Deconvolution") +
  xlab(" ") + ylab("Log(Fold Enrichment)")

pdf(paste0(outputFolder, "targetBased_cellPop_PANTHER_CPM_Deconvolution.pdf"), width = 20, height = 6)
grid.arrange(cpm, deconvolution, nrow = 1)
dev.off()
