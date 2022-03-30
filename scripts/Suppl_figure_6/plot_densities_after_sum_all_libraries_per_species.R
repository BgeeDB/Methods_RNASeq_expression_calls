source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

all_lib_bgee15 <- merge(all_lib_bgee15, sampleInfo, by = "libraryId")
ids <- unique(all_lib_bgee15$speciesId)

## save output
outputFolder <- file.path("./figures/Suppl_figure_6/")
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

listSumSpecies <- list.files(path = "./data/Data_generated_By_Bgee/sum_by_species_15/", pattern = "^sum_abundance_gene_level\\+fpkm\\+intergenic\\+classification_")

for (i in ids) {
  
  summed <- read.table(paste0("./data/Data_generated_By_Bgee/sum_by_species_15/sum_abundance_gene_level+fpkm+intergenic+classification_",i,".tsv"), header=TRUE, sep="\t")
  speciesName <- as.character(unique(all_lib_bgee15$organism[all_lib_bgee15$speciesId == i]))
  numberLibs <- length(all_lib_bgee15$libraryId[all_lib_bgee15$speciesId == i])
  
  pdf(file = paste0(outputFolder, "/distribution_TPM_genic_intergenic_sum_", i, ".pdf"), width = 6, height = 5)
  ## density of log2(TPM) of summed data
  dens <- density(log2(na.omit(summed$tpm) + 10^-6))
  ## genic regions
  dens_genic <- density(c(rep(-30, times=sum(summed$type != "genic")), log2(summed$tpm[summed$type == "genic"] + 10^-6)))
  ## protein-coding genes only
  dens_coding <- density(c(rep(-30, times=sum(!summed$biotype %in% "protein_coding")), log2(summed$tpm[summed$biotype %in% "protein_coding"] + 10^-6)))
  ## intergenic
  dens_intergenic <- density(c(rep(-30, times=sum(summed$type != "intergenic")), log2(summed$tpm[summed$type == "intergenic"] + 10^-6)))
  ## Plot whole distribution
  plot(dens, ylim=c(0, max(c(dens$y, dens_genic$y[dens_genic$x > -15], dens_coding$y[dens_coding$x > -15], dens_intergenic$y[dens_intergenic$x > -15]))*1.1), xlim=c(-23, 21), lwd=2, main=paste0(speciesName, " (", numberLibs, " libraries)"), bty="n", axes=T, xlab="log2(TPM + 10^-6)")
  ## Add subgroups distributions (genic, intergenic, etc):
  ## genic
  lines(dens_genic, col="firebrick3", lwd=2)
  ## protein-coding genes
  lines(dens_coding, col="firebrick3", lwd=2, lty=2)
  ## intergenic
  lines(dens_intergenic, col="dodgerblue3", lwd=2)
  ## legend
  legend("topleft", c(paste0("all (", length(summed[,1]),")"), paste0("genic (", sum(summed$type == "genic"), ")"), paste0("coding (", sum(summed$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed$type == "intergenic"), ")")), lwd=2, col=c("black", "firebrick3", "firebrick3", "dodgerblue3"), lty=c(1, 1, 2, 1), bty="n")
  dev.off()
}

