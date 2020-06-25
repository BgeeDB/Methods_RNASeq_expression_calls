## figure referent to the distribution of calls of expression for Macaca mulatta
## libraries used
library(ggplot2)

pathFiles <- "Files_to_generate_figures/"
outFiles <- "Figure_BgeeCallpaper/"

## collect all data from Bgee 14 and get info about number of libraries
fileData_Bgee14 <- load(paste0(pathFiles, "distribution_comparision_9544.rda"))
fileData_Bgee14 <- all_distribution
librarySize_Bgee14 <- as.character(nrow(dplyr::filter(fileData_Bgee14, approaches == "Intergenic TPM of 1%"))) 
fileData_Bgee14$coding <- as.numeric(fileData_Bgee14$coding)

## collect info from Bgee 14.1
fileData_Bgee14.1 <- read.table(paste0(pathFiles, "Macaca_Bgee14_1.tsv"), header=TRUE, sep="\t")
colnames(fileData_Bgee14.1)[1] <- "libraryId"
subset_Bgee14.1 <- data.frame(fileData_Bgee14.1$proportionCodingPresent, fileData_Bgee14.1$proportionGenicPresent, fileData_Bgee14.1$proportionIntergenicPresent)
colnames(subset_Bgee14.1) <- c("coding", "genic","intergenic")
subset_Bgee14.1$approaches <- "Bgee_14.1"
subset_Bgee14.1$libraryId <- fileData_Bgee14.1$libraryId
librarySize_Bgee14.1 <- (nrow(subset_Bgee14.1)) 

## collect info without Ns
fileData_without_N <- read.table(paste0(pathFiles, "presence_absence_all_samples_without_N.txt"), header=TRUE, sep="\t", comment.char="")
colnames(fileData_without_N)[1] <- "libraryId"
fileData_without_N <- dplyr::filter(fileData_without_N, organism == "Macaca mulatta")
subset_without_N <- data.frame(fileData_without_N$proportionCodingPresent, fileData_without_N$proportionGenicPresent, fileData_without_N$proportionIntergenicPresent)
colnames(subset_without_N) <- c("coding", "genic","intergenic")
subset_without_N$approaches <- "filtered intergenic with 5%"
subset_without_N$libraryId <- fileData_without_N$libraryId
librarySize_without_N <- (nrow(subset_without_N)) 

## collect info without deconvolution
fileData_without_decv <- read.table(paste0(pathFiles, "presence_absence_all_samples_without_decov.txt"), header=TRUE, sep="\t", comment.char="")
colnames(fileData_without_decv)[1] <- "libraryId"
fileData_without_decv <- dplyr::filter(fileData_without_decv, organism == "Macaca mulatta")
subset_without_decv <- data.frame(fileData_without_decv$proportionCodingPresent, fileData_without_decv$proportionGenicPresent, fileData_without_decv$proportionIntergenicPresent)
colnames(subset_without_decv) <- c("coding", "genic","intergenic")
subset_without_decv$approaches <- "No deconvolution"
subset_without_decv$libraryId <- fileData_without_decv$libraryId
librarySize_without_decv <- (nrow(subset_without_decv))

######
## final table to ggplot
allInfo <- rbind(fileData_Bgee14,subset_Bgee14.1, subset_without_N, subset_without_decv)

## plot 
macacaMulatta <- ggplot(allInfo, aes(x=approaches, y=coding)) + 
  geom_boxplot(notch=TRUE) + ylim(0,100) + ylab("Distribution of coding") + xlab(" ") +
  annotate("text", x=c(1,2,3,4,5,6,7), y=100, label= c(librarySize_Bgee14.1, librarySize_without_N, librarySize_Bgee14,librarySize_Bgee14,librarySize_Bgee14,librarySize_Bgee14,librarySize_without_decv )) + 
  ggtitle("Distribution of calls of expressed coding genes")
pdf(paste0(outFiles,"Distribution_Coding_Macaca_mulatta.pdf"), width = 12, height = 5)
print(macacaMulatta)
dev.off()



