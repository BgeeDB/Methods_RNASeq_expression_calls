## Generate the figures for GTeX data (over tissues) and different cutoff (over species)
##  libraries used
library(ggplot2)

pathToFIles <- "Files_to_generate_figures/"
outFiles <- "Figure_BgeeCallpaper/"

Threshold_2 <- read.table(paste0(pathToFIles, "/", "presence_absence_all_samples_TPM-treshold.txt"), header=TRUE, sep="\t", comment.char="")
colnames(Threshold_2)[1] <- "libraryId" 
bgee_withoutDecov <- read.table(paste0(pathToFIles, "/", "presence_absence_all_samples_without_decov.txt"), header=TRUE, sep="\t", comment.char="")
colnames(bgee_withoutDecov)[1] <- "libraryId" 
bgeeThreshold <- read.table(paste0(pathToFIles, "/", "presence_absence_all_samples_14.txt"), header=TRUE, sep="\t", comment.char="")
colnames(bgeeThreshold)[1] <- "libraryId"
bgeeThreshold_without_N <- read.table(paste0(pathToFIles, "/", "presence_absence_all_samples_without_N.txt"), header=TRUE, sep="\t", comment.char="")
colnames(bgeeThreshold_without_N)[1] <- "libraryId"
infoFile <- read.table(paste0(pathToFIles, "/", "info_rnaSeq.tsv"), header=TRUE, sep="\t", comment.char="")
colnames(infoFile)[1] <- "libraryId"

## Threshold = 2
gg_Threshold_2 <- ggplot(data = Threshold_2, mapping = aes(x = organism, y = proportionCodingPresent, label=organism)) +
  geom_boxplot(alpha = 1) +
  labs(x = "", y = "% Protein Coding Genes") +
  scale_y_continuous(limits=c(0,100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file = paste0(outFiles, "/TPM_threshold_allSpecies.pdf"), width = 14, height = 7)
gg_Threshold_2
dev.off()

## Bgee without deconvolution
gg_bgeeWithoutDec <- ggplot(data = bgee_withoutDecov, mapping = aes(x = organism, y = proportionCodingPresent, label=organism)) +
  geom_boxplot(alpha = 1) +
  labs(x = "", y = "% Protein Coding Genes") +
  scale_y_continuous(limits=c(0,100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file = paste0(outFiles, "/Bgee_without_dev_allSpecies.pdf"), width = 14, height = 7)
gg_bgeeWithoutDec
dev.off()

## Bgee Intergenic cut-off
gg_bgee <- ggplot(data = bgeeThreshold, mapping = aes(x = organism, y = proportionCodingPresent, label=organism)) +
  geom_boxplot(alpha = 1) +
  labs(x = "", y = "% Protein Coding Genes") +
  scale_y_continuous(limits=c(0,100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file = paste0(outFiles, "/Bgee_threshold_allSpecies.pdf"), width = 14, height = 7)
gg_bgee
dev.off()

## Bgee Intergenic cut-off without N
gg_bgee_wo_N <- ggplot(data = bgeeThreshold_without_N, mapping = aes(x = organism, y = proportionCodingPresent, label=organism)) +
  geom_boxplot(alpha = 1) +
  labs(x = "", y = "% Protein Coding Genes") +
  scale_y_continuous(limits=c(0,100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file = paste0(outFiles, "/Bgee_threshold_without_N_allSpecies.pdf"), width = 14, height = 7)
gg_bgee_wo_N
dev.off()

## Collect just GTEx data
gtexData <- dplyr::filter(infoFile, experimentId == "SRP012682")

## Gtex with Threshold = 2
gtexThreshold <- merge(Threshold_2, gtexData, by="libraryId")

gtexThreshold_cutoff <- gtexThreshold %>%
  mutate( type=ifelse(uberonName=="blood","Highlighted","Normal")) %>%
  ggplot(aes(x = uberonName, y = proportionCodingPresent, label=uberonName, fill=type, alpha=type)) +
  geom_boxplot(alpha = 1) +
  labs(x = "", y = "% Protein Coding Genes") +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual(values=c("red", "white")) +
  scale_alpha_manual(values=c(1,0.1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "none")

pdf(file = paste0(outFiles, "/GTEx_threshold.pdf"), width = 16, height = 7)
gtexThreshold_cutoff
dev.off()


## Gtex with Bgee cut-off
gtexBgee <- merge(bgeeThreshold, gtexData, by="libraryId")

gtexBgee_cutoff <- gtexBgee %>%
  mutate( type=ifelse(uberonName=="blood","Highlighted","Normal")) %>%
  ggplot(aes(x = uberonName, y = proportionCodingPresent, label=uberonName, fill=type, alpha=type)) +
  geom_boxplot(alpha = 1) +
  labs(x = "", y = "% Protein Coding Genes") +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual(values=c("red", "white")) +
  scale_alpha_manual(values=c(1,0.1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "none")

pdf(file = paste0(outFiles, "/GTEx_Bgee_Cutoff.pdf"), width = 16, height = 7)
gtexBgee_cutoff
dev.off()
