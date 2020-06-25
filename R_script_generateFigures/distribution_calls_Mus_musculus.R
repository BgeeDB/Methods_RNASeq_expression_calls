## figure referent to the distribution of calls of expression for Mus musculus
## libraries used
library(varhandle)
library(ggplot2)

readFile <- load("/Users/sfonseca1/Documents/GitHub_projects/RNASeq_calls/Files_to_generate_figures/distribution_comparision_10090.RDa")
outFiles <- "/Users/sfonseca1/Documents/GitHub_projects/RNASeq_calls/Figure_BgeeCallpaper/"

## collect all data and get info about number of libraries
fileData <- all_distribution
librarySize <- as.character(nrow(dplyr::filter(fileData, approaches == "Intergenic TPM of 1%"))) 
fileData$coding <- unfactor(fileData$coding)

## plot 
musMusculus <- ggplot(fileData, aes(x=approaches, y=coding)) + 
  geom_boxplot(notch=TRUE) + ylim(0,100) + ylab("Distribution of coding") + xlab(" ") +
  annotate("text", x=c(1,2,3,4,5), y=100, label= librarySize) + 
  ggtitle("Distribution of calls of expressed coding genes")
pdf(paste0(outFiles,"Distribution_Coding_Mus_musculus.pdf"), width = 12, height = 5)
print(musMusculus)
dev.off()
