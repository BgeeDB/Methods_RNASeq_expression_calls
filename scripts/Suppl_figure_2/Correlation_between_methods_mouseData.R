source("./scripts//checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- "./figures/Suppl_figure_2/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

outStats <- "./stats_info/Mouse_data/"
if (!dir.exists(outStats)){
  dir.create(outStats)
} else {
  print("Already exists!")
}

####################### collect summary stats to mouse data #####################################3
mouse_0.05_pvalue_decov <- dplyr::filter(allInfo, speciesId == 10090 & approach == "Deconvolution p-value approach" & cutoff == "p-value <= 0.05")
s1 <- summary(mouse_0.05_pvalue_decov$proportionCodingPresent)
mouse_0.05_pvalue_decov <- mouse_0.05_pvalue_decov %>% dplyr::select(libraryId, proportionCodingPresent)
colnames(mouse_0.05_pvalue_decov)[2] <- "pvalue_deconv_0.05"

mouse_0.01_pvalue_decov <- dplyr::filter(allInfo, speciesId == 10090 & approach == "Deconvolution p-value approach" & cutoff == "p-value <= 0.01")
s2 <- summary(mouse_0.01_pvalue_decov$proportionCodingPresent)
mouse_0.01_pvalue_decov <- mouse_0.01_pvalue_decov %>% dplyr::select(libraryId, proportionCodingPresent)
colnames(mouse_0.01_pvalue_decov)[2] <- "pvalue_deconv_0.01"

mouse_0.001_pvalue_decov <- dplyr::filter(allInfo, speciesId == 10090 & approach == "Deconvolution p-value approach" & cutoff == "p-value <= 0.001")
s3 <- summary(mouse_0.001_pvalue_decov$proportionCodingPresent)
mouse_0.001_pvalue_decov <- mouse_0.001_pvalue_decov %>% dplyr::select(libraryId, proportionCodingPresent)
colnames(mouse_0.001_pvalue_decov)[2] <- "pvalue_deconv_0.001"

mouse_Without_decov <- dplyr::filter(allInfo, speciesId == 10090 & approach == "Without deconvolution")
s4 <- summary(mouse_Without_decov$proportionCodingPresent)
mouse_Without_decov <- mouse_Without_decov %>% dplyr::select(libraryId, proportionCodingPresent)
colnames(mouse_Without_decov)[2] <- "pvalue_Withoutdeconv_0.05"

mouse_tpm <- dplyr::filter(allInfo, speciesId == 10090 & approach == "TPM threshold")
s5 <- summary(mouse_tpm$proportionCodingPresent)
mouse_tpm <- mouse_tpm %>% dplyr::select(libraryId, proportionCodingPresent)
colnames(mouse_tpm)[2] <- "tpm"

collectAllStats <- list(s1,s2,s3,s4,s5)
names(collectAllStats) <- c("Decov_0.05", "Decov_0.01", "Decov_0.001", "Without_Decov", "TPM")
capture.output(collectAllStats, file = file.path(outStats,"Statistics_all_approaches_Mouse.tsv"))


### correlation between approaches
get1 <- merge(mouse_0.05_pvalue_decov, mouse_0.01_pvalue_decov, by = "libraryId")
get2 <- merge(get1, mouse_0.001_pvalue_decov, by = "libraryId")
get3 <- merge(get2, mouse_Without_decov, by = "libraryId")
getAll <- merge(get3, mouse_tpm, by = "libraryId")

pdf(file = file.path(outputFolder, paste0("Mouse_Pearson_correlations.pdf")),width = 12, height = 6) 
par(mfrow=c(2,5))
plot(getAll$pvalue_deconv_0.05, getAll$pvalue_deconv_0.001, pch=20,xlab="Deconv p-value <= 0.05", ylab="Deconv p-value <= 0.001")
mtext(paste0("Correlation = " ,round(cor(getAll$pvalue_deconv_0.05, getAll$pvalue_deconv_0.001, method = c("pearson")),2)))

plot(getAll$pvalue_deconv_0.05, getAll$pvalue_deconv_0.01, pch=20,xlab="Deconv p-value <= 0.05", ylab="Deconv p-value <= 0.01")
mtext(paste0("Correlation = " ,round(cor(getAll$pvalue_deconv_0.05, getAll$pvalue_deconv_0.01, method = c("pearson")),2)))

plot(getAll$pvalue_deconv_0.001, getAll$pvalue_deconv_0.01, pch=20,xlab="Deconv p-value <= 0.001", ylab="Deconv p-value <= 0.01")
mtext(paste0("Correlation = " ,round(cor(getAll$pvalue_deconv_0.001, getAll$pvalue_deconv_0.01, method = c("pearson")),2)))

plot(getAll$pvalue_deconv_0.001, getAll$pvalue_Withoutdeconv_0.05, pch=20,xlab="Deconv p-value <= 0.001", ylab="Without deconv.")
mtext(paste0("Correlation = " ,round(cor(getAll$pvalue_deconv_0.001, getAll$pvalue_Withoutdeconv_0.05, method = c("pearson")),2)))

plot(getAll$pvalue_deconv_0.01, getAll$pvalue_Withoutdeconv_0.05, pch=20,xlab="Deconv p-value <= 0.01", ylab="Without deconv.")
mtext(paste0("Correlation = " ,round(cor(getAll$pvalue_deconv_0.01, getAll$pvalue_Withoutdeconv_0.05, method = c("pearson")),2)))

plot(getAll$pvalue_deconv_0.05, getAll$pvalue_Withoutdeconv_0.05, pch=20,xlab="Deconv p-value <= 0.05", ylab="Without deconv.")
mtext(paste0("Correlation = " ,round(cor(getAll$pvalue_deconv_0.05, getAll$pvalue_Withoutdeconv_0.05, method = c("pearson")),2)))

plot(getAll$pvalue_deconv_0.001, getAll$tpm, pch=20,xlab="Deconv p-value <= 0.001", ylab="TPM threshold")
mtext(paste0("Correlation = " ,round(cor(getAll$pvalue_deconv_0.001, getAll$tpm, method = c("pearson")),2)))

plot(getAll$pvalue_deconv_0.01, getAll$tpm, pch=20,xlab="Deconv p-value <= 0.01", ylab="TPM threshold")
mtext(paste0("Correlation = " ,round(cor(getAll$pvalue_deconv_0.01, getAll$tpm, method = c("pearson")),2)))

plot(getAll$pvalue_deconv_0.05, getAll$tpm, pch=20,xlab="Deconv p-value <= 0.05", ylab="TPM threshold")
mtext(paste0("Correlation = " ,round(cor(getAll$pvalue_deconv_0.05, getAll$tpm, method = c("pearson")),2)))

plot(getAll$pvalue_Withoutdeconv_0.05, getAll$tpm, pch=20,xlab="Without deconv.", ylab="TPM threshold")
mtext(paste0("Correlation = " ,round(cor(getAll$pvalue_Withoutdeconv_0.05, getAll$tpm, method = c("pearson")),2)))
dev.off()



