source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- "./figures/Figure_4/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

################################### results TPM >= 2
tpm_threshold <- dplyr::filter(tpm_threshold, libraryId %in% librariesUsed)

##################### results deconvolution p-value <= 0.05
## collect lib + organism info to provide information
libInfo <- tpm_threshold %>% dplyr::select(libraryId, organism)
deconv_pValue_0.05 <- dplyr::filter(deconv_pValue, cutoff == 0.05 & libraryId %in% librariesUsed)
deconv_pValue_0.05 <- merge(deconv_pValue_0.05, libInfo, by = "libraryId")

deconv_pValue_0.05_graph <- ggplot(data = deconv_pValue_0.05, aes(reorder(organism, proportionCodingPresent, median),proportionCodingPresent )) +
  geom_boxplot(alpha = 1) +
  labs(x = "", y = "% Protein Coding Genes") +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle("Deconvolution method")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "none")
pdf(file = paste0(outputFolder, "/Bgee_Deconvolution_allSpecies.pdf"), width = 20, height = 6)
deconv_pValue_0.05_graph
dev.off()


orderByMedian <- aggregate(deconv_pValue_0.05$proportionCodingPresent, list(deconv_pValue_0.05$organism), FUN=median)
orderByMedian <- orderByMedian[order(orderByMedian$x),] 
colnames(orderByMedian) <- c("organism", "ProportionCodingPresent_BgeeMethod")

tpm_threshold <-  merge(tpm_threshold, orderByMedian, by ="organism")


gg_Threshold_2 <- ggplot(tpm_threshold, aes(reorder(organism, ProportionCodingPresent_BgeeMethod), proportionCodingPresent)) +
  geom_boxplot(alpha = 1) +
  labs(x = "", y = "% Protein Coding Genes") +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle("TPM threshold")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "none")
pdf(file = paste0(outputFolder, "/Bgee_TPM_threshold_allSpecies.pdf"), width = 20, height = 6)
gg_Threshold_2
dev.off()

############################ results without deconvolution
without_deconv <- dplyr::filter(without_deconv, libraryId %in% librariesUsed)
without_deconv <- merge(without_deconv, libInfo, by = "libraryId")
without_deconv <-  merge(without_deconv, orderByMedian, by ="organism")

without_deconv_graph <- ggplot(data = without_deconv, aes(reorder(organism, ProportionCodingPresent_BgeeMethod), proportionCodingPresent)) +
  geom_boxplot(alpha = 1) +
  labs(x = "", y = "% Protein Coding Genes") +
  scale_y_continuous(limits=c(0,100)) +
  ggtitle("Without deconvolution")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "none")
pdf(file = paste0(outputFolder, "/Bgee_without_deconv_allSpecies.pdf"), width = 20, height = 6)
without_deconv_graph
dev.off()
