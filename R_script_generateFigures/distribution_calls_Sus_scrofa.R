## figure referent to the distribution of calls for Sus scrofa
## libraries used
library(ggplot2)

pathFiles <- "Files_to_generate_figures/"
outFiles <- "Figure_BgeeCallpaper/"

## collect all data and get info about number of libraries
fileData <- load(paste0(pathFiles, "distribution_comparision_9823.rda"))
fileData <- all_distribution
fileData <- transform(fileData, coding = as.numeric(coding), genic = as.numeric(genic), intergenic = as.numeric(intergenic))
librarySize <- as.character(nrow(dplyr::filter(fileData, approaches == "Intergenic TPM of 1%"))) 

## collect info without deconvolution
fileData_without_decv <- read.table(paste0(pathFiles, "presence_absence_all_samples_without_decov.txt"), header=TRUE, sep="\t", comment.char="")
colnames(fileData_without_decv)[1] <- "libraryId"
fileData_without_decv <- dplyr::filter(fileData_without_decv, organism == "Sus scrofa")
subset_without_decv <- data.frame(fileData_without_decv$proportionCodingPresent, fileData_without_decv$proportionGenicPresent, fileData_without_decv$proportionIntergenicPresent)
colnames(subset_without_decv) <- c("coding", "genic","intergenic")
subset_without_decv$approaches <- "No deconvolution"
subset_without_decv$libraryId <- fileData_without_decv$libraryId

### all info to ggplot
finalInfo <- rbind(fileData, subset_without_decv)

sucScrofa_coding <- ggplot(finalInfo, aes(x=approaches, y=coding)) + 
  geom_boxplot(notch=FALSE) + ylim(0,100) + ylab("Distribution of coding") + xlab(" ") +
  annotate("text", x=c(1,2,3,4,5), y=100, label= librarySize) + 
  ggtitle("Distribution of calls of expressed coding genes")

sucScrofa_genic <- ggplot(finalInfo, aes(x=approaches, y=genic)) + 
  geom_boxplot(notch=FALSE) + ylim(0,100) + ylab("Distribution of genic") + xlab(" ") +
  annotate("text", x=c(1,2,3,4,5), y=100, label= librarySize) + 
  ggtitle("Distribution of genic expressed")

sucScrofa_intergenic <- ggplot(finalInfo, aes(x=approaches, y=intergenic)) + 
  geom_boxplot(notch=FALSE) + ylim(0,100) + ylab("Distribution of intergenic") + xlab(" ") +
  annotate("text", x=c(1,2,3,4,5), y=100, label= librarySize) + 
  ggtitle("Distribution of intergenic expressed")

pdf(paste0(outFiles,"All_Distributions_Sus_scrofa.pdf"), width = 12, height = 5)
print(sucScrofa_coding)
print(sucScrofa_genic)
print(sucScrofa_intergenic)
dev.off()


