source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- "./figures/Figure_7/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

## reference intergenic regions
refInt_mouse <- read.table("./data/Reference_intergenic_regions/10090_coordinates.tsv", header=TRUE, sep="\t")

##### one cell
pickOneCell <- oneCell %>% dplyr::select(gene_id, AAACCTGAGCAGCCTC,biotype, type, cellTypeName, cellTypeId)

dens <- density(log2(pickOneCell$AAACCTGAGCAGCCTC+1e-6), na.rm=T)
dens_genic <- density(log2(pickOneCell$AAACCTGAGCAGCCTC[pickOneCell$type == "genic"]+1e-6), na.rm=T)
dens_genic$y <- dens_genic$y * nrow(dplyr::filter(pickOneCell, type == "genic")) / length(pickOneCell$AAACCTGAGCAGCCTC)
dens_coding <- density(log2(pickOneCell$AAACCTGAGCAGCCTC[pickOneCell$biotype == "protein_coding"]+1e-6), na.rm=T)
dens_coding$y <- dens_coding$y * nrow(dplyr::filter(pickOneCell, biotype == "protein_coding")) / length(pickOneCell$AAACCTGAGCAGCCTC)
dens_intergenic <- density(log2(pickOneCell$AAACCTGAGCAGCCTC[pickOneCell$type == "intergenic"]+1e-6), na.rm=T)
dens_intergenic$y <- dens_intergenic$y * nrow(dplyr::filter(pickOneCell, type == "intergenic")) / length(pickOneCell$AAACCTGAGCAGCCTC)

## ref intergenic coordinates
extract_Refintergenic <- dplyr::filter(pickOneCell, gene_id %in% refInt_mouse$chr_start_end)
dens_Refintergenic <- density(log2(extract_Refintergenic$AAACCTGAGCAGCCTC+1e-6), na.rm=T)
dens_Refintergenic$y <- dens_Refintergenic$y * nrow(extract_Refintergenic) / length(pickOneCell$AAACCTGAGCAGCCTC)

pdf(paste0(outputFolder, "targetBased_OneCell.pdf"), width = 6, height = 6)
## plot density one cell
plot(dens, lwd=2, main="One B-cell", xlab="Log2(CPM)")
mtext(" Cell id: AAACCTGAGCAGCCTC")
lines(dens_genic,col="red", lwd=2)
lines(dens_coding,col="indianred", lwd=2)
lines(dens_intergenic,col="darkblue", lwd=2)
lines(dens_Refintergenic,col="cyan", lwd=2)
legend("topright", legend = c("All", "Genic" ,"Protein coding", "Intergenic", "Ref_Intergenic"), col=c("black", "red", "indianred", "darkblue", "cyan"),lty=c(1,1,1,1,1), lwd=2, bty = "n")
dev.off()

## zoom on the right part of the plot
pdf(paste0(outputFolder, "targetBased_OneCell_ZOOM.pdf"), width = 6, height = 6)
## plot density one cell
plot(dens, lwd=2, main="One B-cell", xlab="Log2(CPM)", xlim=c(5,17), ylim=c(0,0.015))
mtext(" Cell id: AAACCTGAGCAGCCTC")
lines(dens_genic,col="red", lwd=2)
lines(dens_coding,col="indianred", lwd=2)
lines(dens_intergenic,col="darkblue", lwd=2)
lines(dens_Refintergenic,col="cyan", lwd=2)
legend("topright", legend = c("All", "Genic" ,"Protein coding", "Intergenic", "Ref_Intergenic"), col=c("black", "red", "indianred", "darkblue", "cyan"),lty=c(1,1,1,1,1), lwd=2, bty = "n")
dev.off()


#### population of B-cells
pdf(paste0(outputFolder, "targetBased_BCellPop.pdf"), width = 6, height = 6)
dens <- density(log2(bCellPop$CPM+1e-6), na.rm=T)
dens_genic <- density(log2(bCellPop$CPM[bCellPop$type == "genic"]+1e-6), na.rm=T)
dens_genic$y <- dens_genic$y * nrow(dplyr::filter(bCellPop, type == "genic")) / length(bCellPop$CPM)
dens_coding <- density(log2(bCellPop$CPM[bCellPop$biotype == "protein_coding"]+1e-6), na.rm=T)
dens_coding$y <- dens_coding$y * nrow(dplyr::filter(bCellPop, biotype == "protein_coding")) / length(bCellPop$CPM)
dens_intergenic <- density(log2(bCellPop$CPM[bCellPop$type == "intergenic"]+1e-6), na.rm=T)
dens_intergenic$y <- dens_intergenic$y * nrow(dplyr::filter(bCellPop, type == "intergenic")) / length(bCellPop$CPM)

## ref intergenic coordinates
extract_Refintergenic <- dplyr::filter(bCellPop, gene_id %in% refInt_mouse$chr_start_end)
dens_Refintergenic <- density(log2(extract_Refintergenic$CPM+1e-6), na.rm=T)
dens_Refintergenic$y <- dens_Refintergenic$y * nrow(extract_Refintergenic) / length(bCellPop$CPM)

## plot density population of B-cells
plot(dens, lwd=2, main="B-Cell population", xlab="Log2(CPM)")
mtext(" SRX6060813 - CL-0000236")
lines(dens_genic,col="red", lwd=2)
lines(dens_coding,col="indianred", lwd=2)
lines(dens_intergenic,col="darkblue", lwd=2)
lines(dens_Refintergenic,col="cyan", lwd=2)
legend("topright", legend = c("All", "Genic" ,"Protein coding", "Intergenic", "Ref_Intergenic"), col=c("black", "red", "indianred", "darkblue", "cyan"),lty=c(1,1,1,1,1), lwd=2, bty = "n")
dev.off()




