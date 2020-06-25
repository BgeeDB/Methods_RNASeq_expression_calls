## generate spearman correlation graphic for gtex data

## libraries used
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

pathFiles <- "Files_to_generate_figures/"
outFiles <- "Figure_BgeeCallpaper/"

infoRNASeq <- fread(paste0(pathFiles, "info_rnaSeq.tsv"))
colnames(infoRNASeq)[1] <- "libraryId"

deconvolutionFile <- read.table(paste0(pathFiles, "presence_absence_all_samples_14.txt"), header=TRUE, sep="\t", comment.char="")
colnames(deconvolutionFile)[1] <- "libraryId"
tpmThreshold <-  read.table(paste0(pathFiles, "presence_absence_all_samples_TPM-treshold.txt"), header=TRUE, sep="\t", , comment.char="")
colnames(tpmThreshold)[1] <- "libraryId"

## extract GTeX experiment
gtex <- "SRP012682"
selectLibrariesExperiment <- dplyr::filter(infoRNASeq, experimentId == gtex)
selectUberonName <- data.frame(selectLibrariesExperiment$libraryId, selectLibrariesExperiment$uberonName)
colnames(selectUberonName) <- c("libraryId", "Uberon_Name")
selectLibraries <- unlist(selectLibrariesExperiment[,1], use.names=FALSE) 

deconvolutionFile <- filter(deconvolutionFile, libraryId %in% selectLibraries)
deconvolutionFile <- merge(deconvolutionFile, selectUberonName, by = "libraryId")

tpmThreshold <- filter(tpmThreshold, libraryId %in% selectLibraries)
tpmThreshold <- merge(tpmThreshold, selectUberonName, by = "libraryId")

### plot data
finalTable <- data.frame(deconvolutionFile$libraryId, deconvolutionFile$proportionCodingPresent, tpmThreshold$proportionCodingPresent, tpmThreshold$Uberon_Name)
colnames(finalTable) <- c("libraryId", "PC_deconvolution", "PC_TPM", "Uberon_Name")

spearmanPlot <- ggplot(finalTable,
            aes(x=PC_deconvolution, y=PC_TPM)) +
  geom_point(color = dplyr::case_when(finalTable$Uberon_Name == "blood" ~ "indianred",TRUE ~ "gray"), size = 2) +
  ylim(0,100) +xlim(0,100) + ggtitle("GTeX libraries")+
  ylab(expression(paste("% coding present TPM >= 2"))) +
  xlab("% coding present deconvolution") +
  stat_cor(method = "pearson", label.x = 3, label.y = 95) 
pdf(paste0(outFiles, "Spearman_correlation_GTeX.pdf"), width = 8, height = 6)
print(spearmanPlot)
dev.off()
