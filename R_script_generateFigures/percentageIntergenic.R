# Code to generate the figure percentage of intergenic per species
##  libraries used
library(ggplot2)

pathToFIles <- "Files_to_generate_figures/"
outFiles <- "Figure_BgeeCallpaper/"

readFile <- read.table(paste0(pathToFIles, "percentage_intergenic_species.tsv"), header=TRUE, sep="\t", comment.char="")
readFile <- readFile[order(readFile$species),]

percentageIntergenic <- ggplot(data=readFile, aes(x=species, y=Percentage_Ns, fill=Region)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme(axis.text.x = element_text(angle = 90)) +
 scale_fill_brewer(palette="Paired") + ylab("% of Ns") + xlab(" ") + ylim(0,100) +
  ggtitle("Percentage of intergenic per species")
pdf(paste0(outFiles,"Percentage_Intergenic.pdf"), width = 12, height = 5)
print(percentageIntergenic)
dev.off()
