## riboSeq genes Ids identified as TP and TN for Mouse liver data
TP_TN_file <- "./data/benchmark_files/TP_TN_usingJust_ProteinCoding_TP_filter_TN_others.tsv"

## Epigenetic markers for Human lung data
TP_TN_human_lung <- "./data/benchmark_files/epigeneticMarkers.tsv"

## Testis markers from developmental studies for Drosophila data
TP_TN_drosophila_testis <- "./data/benchmark_files/markers_fromZigzag_Drosophila.tsv"


## Input folders for benchmark Mouse liver
folderFiles_AllIntergenic <- "./data/MouseData_allLibraries_for benchmark/just_Mouse_liver_samples_with_all_intergenic/"
folderFiles_RefIntergenic <- "./data/MouseData_allLibraries_for benchmark/just_Mouse_liver_RefIntergenic/"

## Input folders for benchmark Human lung
folderHumanLungAllInt <- "./data/HumanData_allLibraries_for_benchmark/Human_lung_AllIntergenic/"
folderHumanLungRefInt <- "./data/HumanData_allLibraries_for_benchmark/Human_lung_RefIntergenic/"

## Input folders for benchmark Drosophila testis
folderDrosTestisAllInt <- "./data/DrosophilaData_allLibraries_for_benchmark/Drosophila_testis_AllIntergenic/"
folderDrosTestisRefInt <- "./data/DrosophilaData_allLibraries_for_benchmark/Drosophila_testis_RefIntergenic/"


## TPM input files (collected from Bgee 15) used in zigzag method
humanData <- read.table("./complementary_analysis/zigzag/Input_files_zigzag/Human_all_lung_samples_TPM_Bgee15.tsv", header=TRUE, sep="\t")
flyData <- read.table("./complementary_analysis/zigzag/Input_files_zigzag/Drosophila_all_testis_samples_TPM_Bgee15.tsv", header=TRUE, sep="\t")