source("./scripts/checkLibraries.R")

## Secies name and correspondent species ID
speciesInfoCollected <- read.table("./data/Data_generated_By_Bgee/SpeciesName_SpeciesID_id.tsv", header=TRUE, sep="\t")

## all libraries annotated in Bgee 15.0
all_lib_bgee15 <- read.table("./data/Data_generated_By_Bgee/all_bgee15_RNAseq_lib.txt", header=TRUE, sep="\t")
librariesUsed <- all_lib_bgee15$rnaSeqLibraryId
colnames(all_lib_bgee15) <- "libraryId"

## info about each bulk rna-seq library 
sampleInfo <- read.table("./data/Data_generated_By_Bgee/rna_seq_sample_info.txt", header=TRUE, sep="\t")

## annotation of all RNA-Seq libraries
rnaSeqLib <- read.table("./data/Data_generated_By_Bgee/rna_seq_libraries.tsv", header=TRUE, sep="\t")
colnames(rnaSeqLib)[1] <- "libraryId"

## sum file of all libraries for a target species to select the reference intergenic regions (Mouse and Pike)
sum_by_species_mouse <- "./data/Data_generated_By_Bgee/sum_by_species_15/sum_abundance_gene_level+fpkm+intergenic+classification_10090.tsv"
sum_by_species_pike <- "./data/Data_generated_By_Bgee/sum_by_species_15/sum_abundance_gene_level+fpkm+intergenic+classification_8010.tsv"

## reference intergenic regions from Bgee 15 
pike_refIntergenic <- read.table("./data/Reference_intergenic_regions/8010_coordinates.tsv", header=TRUE, sep="\t")

## calls of expression using different approaches (tmp threshold >= 2, deconvolution using different p-values cut-off or without deconvolution with p-value cutoff <= 0.05)
tpm_threshold <- read.table("./data/Calls_expression/Calls_perLibrary_TPM_2_Genic_ProteinCoding.tsv", header=TRUE, sep="\t")
deconv_pValue <- read.table("./data/Calls_expression/summary_Stats_All_Libraries_different_pValues.tsv", header=TRUE, sep="\t")
without_deconv <- read.table("./data/Calls_expression/summary_Stats_All_Libraries_NO_decov.tsv", header=TRUE, sep="\t")

## proportion of protein coding genes per library using all methods (file contain info above)
allInfo <- read.table("./data/Calls_expression/CodingGenes_callPresent_perLibrary_allApproaches.tsv", header = T, sep="\t")


## summary stats file with proportion of coding genes for bulk RNA-Seq and scRNA-Seq full length data.
RNA_scRNA <- read.table("./data/Calls_expression/RNA+scRNA_Calls_PresenceAbsence_all_samples_Bgee15.txt", header=TRUE, sep="\t")

## calls of expression using different reference intergenic methods (pValue or qValue) to call present and absent genes in (Human lung, Drosophila testis and Mouse liver samples annotated in Bgee 15)
callExpression <- read.table("./data/Calls_expression/summary_Stats_All_Libraries_Used_Combining_process_for_pValue_and_qValue_method.tsv", header=TRUE, sep="\t")

## Droplet-based analysis (one cell + cell population)
oneCell <- fread("./data/Droplet_scRNA-Seq/All_FINAL_RESULTS/SRX6060813/Normalized_Counts_B-cell_CL-0000236.tsv", header=TRUE, sep="\t")
oneCell$gene_id <- gsub("-", "_", oneCell$gene_id)

## B cell population from target-based 
bCellPop <- read.table("./data/Droplet_scRNA-Seq/All_FINAL_RESULTS/SRX6060813/Calls_cellPop_SRX6060813_B-cell_CL-0000236_genic+intergenic.tsv", header=TRUE, sep = "\t")
bCellPop$gene_id <- gsub("-", "_", bCellPop$gene_id)

## info of mouse
mouse_tx2gene <- "./data/Species_information/Mus_musculus.GRCm38.tx2gene"
mouse_gene2biotype <- "./data/Species_information/Mus_musculus.GRCm38.gene2biotype"

## reference intergenic from Bgee 15 for mouse 
refInt_mouse <- read.table("./data/Reference_intergenic_regions/10090_coordinates.tsv", header=TRUE, sep="\t")

## results all Droplet-based samples
targetBased <- fread("./data/Droplet_scRNA-Seq/All_FINAL_RESULTS/All_cellPopulation_stats_10X.tsv")

## scRNA-Seq library file
scRNAlibraryInfo <- read.table("./data/Data_generated_By_Bgee/scRNASeqLibrary.tsv", header=TRUE, sep="\t")
colnames(scRNAlibraryInfo)[1] <- "libraryId"



