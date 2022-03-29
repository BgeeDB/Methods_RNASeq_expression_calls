## Usage:
## R CMD BATCH --no-save --no-restore '--args librariesFile="librariesFile.tsv" species="speciesID" geneLengthFile="geneLengthFile.tsv" outputFolder="output"' zigzag_analysis.R zigzag_analysis.Rout

library(data.table)
library(zigzag)
library(dplyr)

## reading arguments
cmd_args = commandArgs(TRUE);
print(cmd_args)
if( length(cmd_args) == 0 ){ stop("no arguments provided\n") } else {
  for( i in 1:length(cmd_args) ){
    eval(parse(text=cmd_args[i]))
  }
}

## checking if all necessary arguments were passed....
command_arg <- c("librariesFile","species","geneLengthFile","outputFolder")
for( c_arg in command_arg ){
  if( !exists(c_arg) ){
    stop( paste(c_arg,"command line argument not provided\n") )
  }
}

## Read file with TPM per library. If file not exists, script stops
if( file.exists(librariesFile) ){
  libraries <- read.table(librariesFile, h=T, sep="\t")
} else {
  stop( paste("The libraries file not found [", librariesFile, "]\n"))
}

## Read geneLengthFile file. If file not exists, script stops
if( file.exists(geneLengthFile) ){
  geneLength <- read.table(geneLengthFile, header = TRUE, sep="\t")
} else {
  stop( paste("The file not found [", geneLengthFile, "]\n"))
}

################### Create gene-length file per species ###########################
# geneLengthInfo <- function(geneLength, species){
#   
#   cat("Preparing gene length for the species: ", species, "\n")
#   
#   ## select just information referent to the species
#   dataSpecies <- dplyr::filter(geneLength, V1 == species)
#   
#   ## creat a table per gene length (this step is a bit longer!!!)
#   finalTable <- c()
#   
#   for (i in unique(dataSpecies$V3)) {
#     
#     ID <- i
#     dataInfo <- dplyr::filter(dataSpecies, V3 == paste0(ID))
#     collectLengthGene <-  sum(dataInfo$V4) / length(dataInfo$V4)
#     collect <- data.frame(ID, collectLengthGene)
#     finalTable <- rbind(finalTable, collect)
#   }
#   return(finalTable)
# }
# if (file.exists(file.path(outputFolder, paste0("geneLength_",species,".tsv"))) == TRUE){
#   message("Don't need to be generated.", "\n")
#   geneLength_Info <- read.table(file.path(outputFolder, paste0("geneLength_",species,".tsv")), header=TRUE, sep="\t")
# } else {
#   geneLength_Info <- geneLengthInfo(geneLength=geneLength, species=species)
#   write.table(geneLength_Info, file = file.path(outputFolder, paste0("geneLength_",species,".tsv")), row.names = FALSE, quote = FALSE, sep = "\t")
# }

#########################  GTFtool to get gene mean for each ID (alternative to the code above that use Bgee information file)
## use gene mean calculated using the GTFtool for each species
geneLength <- geneLength[c(1,2)]
## order per geneID
geneLength <- geneLength[order(geneLength$gene) , ]
## keep just gene ID that are in the libraries 
geneLength <- geneLength[ geneLength$gene %in% libraries$id, ]
## provide gene ID as row.name
row.names(geneLength) <- geneLength$gene
geneLength$gene <- NULL
geneLength$mean <- as.numeric(geneLength$mean)


################### Start analysis of the data with zigzag ###########################
cat("Start running zigzag method", "\n")

pdf(paste0(outputFolder, "/Distribution_samples_",species,".pdf"))
## plot of the density distribution
plot(density(log(libraries[,2])), main =paste0("Distribuiton of samples - ", species), xlab = "log Expression", ylim = c(0, 0.25))
for(i in 2:ncol(libraries)) lines(density(log(libraries[,i])))
abline(v = c(1,4), lwd = c(2, 1), col = c("darkblue", "indianred"))
dev.off()

## load the data into zigzag object
## with active expression component between boundaries 1 and 4 (default)
ptm <- proc.time()
rownames(libraries) <- libraries$id
libraries$id <- NULL
my_zigzag <- zigzag$new(data = libraries, gene_length = geneLength, 
                        output_directory = paste0(outputFolder,"zigzag_analysis_output"), 
                        num_active_components = 2, threshold_a = c(1,4), threshold_i = 1)
proc.time()-ptm

cat("Burning step...", "\n")
ptm <- proc.time()
## Looking the burnin (convergence). NOTE: sample_frequency and ngen can make the process more slow (but we increse the ngen to 50000)
burnin <- my_zigzag$burnin(sample_frequency = 50, progress_plot = FALSE, write_to_files = TRUE, ngen=50000)
proc.time()-ptm

cat("MCMC chain...", "\n")
## Run the MCMC analysis --> in this step we can decrease the number of generations but this can influence the MCMC chain, so keep as suggested by zigzag. 
## we don't use run_posterior_predictive = TRUE because if large number of samples the process is even more slow.
ptm <- proc.time()
mcmc <- my_zigzag$mcmc(sample_frequency = 50, ngen = 50000, mcmcprefix = paste0(species))
proc.time()-ptm

