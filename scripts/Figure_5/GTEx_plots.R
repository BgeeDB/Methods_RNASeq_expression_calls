source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- "./figures/Figure_5/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

homoRnaSeq <- Bgee$new(species = "Homo_sapiens", dataType = "rna_seq", release = "15.0")
homoAnnot <- getAnnotation(homoRnaSeq)

homoAnnot <- read.table("./Homo_sapiens_Bgee_15_0/Homo_sapiens_RNA-Seq_libraries.tsv", header=TRUE, sep="\t")
homoAnnot <- dplyr::filter(homoAnnot, Experiment.ID == "SRP012682")
homoAnnot <- data.frame(homoAnnot$Experiment.ID ,homoAnnot$Library.ID, homoAnnot$Anatomical.entity.name)
colnames(homoAnnot)[2] <- "libraryId"

libraries_GTEX <- homoAnnot$libraryId

## results deconvolution p-value <= 0.05
deconv_pValue_0.05 <- dplyr::filter(deconv_pValue, cutoff == 0.05 & libraryId %in% libraries_GTEX)
deconv_pValue_0.05 <- merge(deconv_pValue_0.05, homoAnnot, by = "libraryId")
deconv_pValue_0.05 <- deconv_pValue_0.05[order(deconv_pValue_0.05$homoAnnot.Anatomical.entity.name),]


gtex_pValue_BgeeApproach <- deconv_pValue_0.05 %>%
  mutate( type=ifelse(homoAnnot.Anatomical.entity.name=="blood","Highlighted","Normal")) %>%
  ggplot( aes(reorder(homoAnnot.Anatomical.entity.name, proportionCodingPresent, median), proportionCodingPresent, label=homoAnnot.Anatomical.entity.name, fill=type, alpha=type)) +
  geom_boxplot(alpha = 1) +
  labs(x = "", y = "% Protein Coding Genes") +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual(values=c("red", "white")) +
  scale_alpha_manual(values=c(1,0.1)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "none")

pdf(file = paste0(outputFolder, "/GTEx_Deconv_pValue_0.05_BgeeApproach.pdf"), width = 16, height = 7)
gtex_pValue_BgeeApproach
dev.off()

orderByMedian <- aggregate(deconv_pValue_0.05$proportionCodingPresent, list(deconv_pValue_0.05$homoAnnot.Anatomical.entity.name), FUN=median)
orderByMedian <- orderByMedian[order(orderByMedian$x),] 
colnames(orderByMedian) <- c("homoAnnot.Anatomical.entity.name", "ProportionCodingPresent_BgeeMethod")

### results TPM >= 2
tpm_threshold <- dplyr::filter(tpm_threshold, libraryId %in% libraries_GTEX)
tpm_threshold <- merge(tpm_threshold, homoAnnot, by = "libraryId")
tpm_threshold <- tpm_threshold[order(tpm_threshold$homoAnnot.Anatomical.entity.name),]
tpm_threshold <-  merge(tpm_threshold, orderByMedian, by ="homoAnnot.Anatomical.entity.name")

gtexThreshold_cutoff <- tpm_threshold %>%
  mutate( type=ifelse(homoAnnot.Anatomical.entity.name=="blood","Highlighted","Normal")) %>%
  ggplot(aes(reorder(homoAnnot.Anatomical.entity.name, ProportionCodingPresent_BgeeMethod), y = proportionCodingPresent, label=homoAnnot.Anatomical.entity.name, fill=type, alpha=type)) +
  geom_boxplot(alpha = 1) +
  labs(x = "", y = "% Protein Coding Genes") +
  scale_y_continuous(limits=c(0,100)) +
  scale_fill_manual(values=c("red", "white")) +
  scale_alpha_manual(values=c(1,0.1)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "none")

pdf(file = paste0(outputFolder, "/GTEx_TPM_threshold.pdf"), width = 16, height = 7)
gtexThreshold_cutoff
dev.off()

### Pearson correlation between TPM and deconv p-value ≤ 0.05 
tpm_threshold <- tpm_threshold %>% dplyr::select(libraryId, proportionCodingPresent, homoAnnot.Anatomical.entity.name)
colnames(tpm_threshold)[2] <- "proportionCodingPresent_TPM"
deconv_pValue_0.05 <- deconv_pValue_0.05 %>% dplyr::select(libraryId, proportionCodingPresent, homoAnnot.Anatomical.entity.name)
colnames(deconv_pValue_0.05)[2] <- "proportionCodingPresent_pValue"
getAllInfo <- merge(tpm_threshold, deconv_pValue_0.05, by = c("libraryId", "homoAnnot.Anatomical.entity.name"))

pearsonPlot <- ggplot(getAllInfo,
                      aes(x=proportionCodingPresent_pValue, y=proportionCodingPresent_TPM)) +
  geom_point(color = dplyr::case_when(getAllInfo$homoAnnot.Anatomical.entity.name == "blood" ~ "indianred",TRUE ~ "gray"), size = 2) +
  ylim(0,100) +xlim(0,100) + ggtitle("GTEx libraries")+
  ylab(expression(paste("% coding present TPM >= 2"))) +
  xlab("% coding present pvalue method") +
  stat_cor(method = "pearson", label.x = 3, label.y = 95)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "none")

pdf(paste0(outputFolder, "Pearson_correlation_GTEx.pdf"), width = 8, height = 6)
print(pearsonPlot)
dev.off()
