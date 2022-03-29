source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

## save output
outputFolder <- file.path("./figures/Suppl_figure_7/")
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

################################## get correlation from GTEX experiment without blood samples
rnaSeqLibBlood <- rnaSeqLib[rnaSeqLib$experimentId == "SRP012682" & rnaSeqLib$uberonName == "blood",]
bloodLib <- rnaSeqLibBlood$libraryId

## get just gtex samples
rnaSeqLibGtex <- rnaSeqLib$libraryId[rnaSeqLib$experimentId == "SRP012682"]

bgeeApproach <- deconv_pValue[deconv_pValue$cutoff == "0.05",]
bgeeApproach <- dplyr::filter(bgeeApproach, !libraryId %in% bloodLib)
bgeeApproach <- dplyr::filter(bgeeApproach, libraryId %in% rnaSeqLibGtex)

tpm_threshold <- dplyr::filter(tpm_threshold, !libraryId %in% bloodLib)
tpm_threshold <- dplyr::filter(tpm_threshold, libraryId %in% rnaSeqLibGtex)

allInfo_withou_blood_gtex <- merge(bgeeApproach, tpm_threshold, by="libraryId")

pearsonPlotWithoutBlood <- ggplot(allInfo_withou_blood_gtex,
                                  aes(x=proportionCodingPresent.x, y=proportionCodingPresent.y))+
  geom_point()+
  ylim(0,100) +xlim(0,100) + ggtitle("GTEx libraries without blood samples")+
  ylab(expression(paste("% coding present TPM >= 2"))) +
  xlab("% coding present deconvolution") +
  stat_cor(method = "pearson", label.x = 3, label.y = 95)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "none")

pdf(paste0(outputFolder, "Pearson_correlation_GTEx_without_blood_samples.pdf"), width = 8, height = 6)
print(pearsonPlotWithoutBlood)
dev.off()
