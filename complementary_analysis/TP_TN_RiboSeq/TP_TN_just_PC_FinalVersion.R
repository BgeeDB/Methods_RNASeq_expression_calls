source("./scripts/checkLibraries.R")

pathFolder <- "./complementary_analysis/TP_TN_RiboSeq/RibosomeFootprint/"

## select all files from all time points
files <- list.files(pathFolder, pattern = "*.rpkm.txt$", full.names = TRUE)
AllFiles <- lapply(files, read.delim, header=FALSE, sep="\t")
DATA <- do.call("cbind", AllFiles)
select_info <- DATA[,1]
select_RPKM = DATA[, grepl("^V2", names( DATA))]

finalData <- data.frame(select_info, select_RPKM)
colnames(finalData)[1] <- "id"

## get Info gene to biotype
gene2Biotype <- read.table("./data/Species_information/Mus_musculus.GRCm38.gene2biotype",header=TRUE, sep="\t")
allInfo <- merge(finalData, gene2Biotype, by ="id")

selectJust_PC <- dplyr::filter(allInfo, biotype == "protein_coding")
## remove 2 last columns (biotype and type)
selectJust_PC <- selectJust_PC[, -c(26:27)]

### collect the TRUE POSITIVES
## select genes with RPKM > 1 in all samples (to create the TP list)
high_1_allSamples <- selectJust_PC[apply(selectJust_PC[,2:ncol(selectJust_PC)] > 1, 1, all), ]
high_1_allSamples$call <- "present"
high_1_allSamples <- high_1_allSamples %>% dplyr::select(id, call)
## select genes with RPKM > 5 in at least one sample (to create the TP list)
high_5_oneSample <- selectJust_PC[apply(selectJust_PC[,2:ncol(selectJust_PC)] > 5, 1, any), ]
high_5_oneSample$call <- "present"
high_5_oneSample <- high_5_oneSample %>% dplyr::select(id, call)
## Union of the 2 data.frame to create the list of True Positive genes
tp <- union(high_1_allSamples$id, high_5_oneSample$id)
tp <- as.data.frame(tp)
tp$call <- "present"
colnames(tp)[1] <- "id"

## select TN
selectAll_PC <- dplyr::filter(gene2Biotype, biotype == "protein_coding")
tn <- dplyr::filter(selectAll_PC, !(id %in% tp$id) )
tn$call <- "absent"
tn <- tn[,c(1,4)]

allInfo_TP_TN_RiboSeq <- rbind(tp, tn)
write.table(allInfo_TP_TN_RiboSeq, file = "./complementary_analysis/TP_TN_RiboSeq/TP_TN_usingJust_ProteinCoding_TP_filter_TN_others.tsv", sep="\t", quote = FALSE, row.names = FALSE)
