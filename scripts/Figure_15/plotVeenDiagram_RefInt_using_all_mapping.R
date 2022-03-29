source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- file.path("./figures/Figure_15/")
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

### collect the sum files from the different mapping or pseudo-mapping tool
sumFile_salmon <- read.table("./data/Different_alignments_PIKE/pike_folder_samples_Quants_Salmon/sum_abundance_gene_level+fpkm+intergenic+classification_8010_salmon.tsv", header=TRUE, sep="\t")
sumFile_hisat2 <- read.table("./data/Different_alignments_PIKE/pike_folder_samples_Quants_Hisat2/sum_abundance_gene_level+fpkm+intergenic+classification_8010_Hisat2.tsv", header=TRUE, sep="\t")


#extract reference intergenic where TPM < TPM intergenic 2
select_intergenic <- dplyr::filter(sumFile_salmon, type == "intergenic", classification == "intergenic_2")
max_TPM <- max(select_intergenic$tpm)
refIntergenic <- dplyr::filter(sumFile_salmon, type == "intergenic", tpm <= max_TPM)
select_intergenicRef_salmon <- refIntergenic[,1]

#extract reference intergenic from Hisat2 (select gaussian 2)
select_intergenic <- dplyr::filter(sumFile_hisat2, type == "intergenic", classification == "intergenic_2")
max_TPM <- max(select_intergenic$tpm)
refIntergenic <- dplyr::filter(sumFile_hisat2, type == "intergenic", tpm <= max_TPM)
select_intergenicRef_hisat2 <- refIntergenic[,1]


### cross this information with ref intergenic selected with kallisto
pike_refIntergenic <- unlist(pike_refIntergenic$chr_start_end)

pdf(paste0(outputFolder, "/Compare_refInt_Kallisto_salmon_Hisat2_pike.pdf"), width = 6, height = 6)
allData <- list(select_intergenicRef_salmon, select_intergenicRef_hisat2, pike_refIntergenic)
vp_alpha <- venn.diagram(allData, fill = c("red", "cyan", "violet"), alpha = 0.3, filename = NULL, 
                         category.names=c("Salmon", "Hisat2","Kallisto"), cex = 1.5)
grid.draw(vp_alpha)
dev.off()
