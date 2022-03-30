source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

speciesId <- "Esox lucius"
gaussianSelected <- 2

outputFolder <- file.path("./figures/Figure_3/")
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

summed_filtered <- read.table(sum_by_species_pike, header=TRUE, sep="\t")
## get_cutoff
select_intergenic <- dplyr::filter(summed_filtered, type == "intergenic", classification == paste0("intergenic_",gaussianSelected)) 
max_TPM <- max(select_intergenic$tpm)
refIntergenic <- dplyr::filter(summed_filtered, type == "intergenic", tpm <= max_TPM)
select_intergenicRef <- refIntergenic[,1]

pdf(file = file.path(outputFolder, paste0("Gaussians_Density_Distributions_",speciesId,".pdf")), width = 10, height = 6)
par(mfrow=c(1,2))
summed_filtered <- summed_filtered[summed_filtered$tpm > 10^-6, ]
summed_filteredIntergenic <- summed_filtered[summed_filtered$type == "intergenic", ]
set.seed(123)
mc <- Mclust(log2(summed_filteredIntergenic$tpm))
mod2 <- densityMclust(log2(summed_filteredIntergenic$tpm), plot = FALSE)
cols <- brewer.pal(mod2$G,'Set2')
for (i in 1 : mod2$G){
  meanData <- mod2$parameters$mean[i]
  sdData <- sqrt(mod2$parameters$variance$sigmasq[i])
  propData <- mod2$parameters$pro[i]
  d <- density(log2(summed_filteredIntergenic$tpm[summed_filteredIntergenic$type == "intergenic" & summed_filteredIntergenic$classification == paste0("intergenic_",i)]), bw = 1)
  d1 <- dnorm(d$x, meanData, sdData)*propData
  if(i == 1){
    plot(d$x, d1, col=cols[i], lwd = 3, ylim=c(0,0.2), type="l", main=speciesId, xlab="log2(TPM)", ylab="Density", bty="n")
    mtext("Gaussians")
  } else {
    lines(d$x, d1, col=cols[i], lwd = 3, ylim=c(0,0.20), type="l") 
  }
}


## Plot the density of the original data, and the density of regions classified to different gaussians
dens <- density(log2(summed_filtered$tpm))
## Plot whole distribution
dens_coding <- density(log2(summed_filtered$tpm[summed_filtered$biotype %in% "protein_coding"]))
dens_coding$y <- dens_coding$y * sum(summed_filtered$biotype %in% "protein_coding") / length(summed_filtered$tpm)
plot(dens_coding, col="firebrick3", lwd=2, lty=2, ylim=c(0,0.15), main=speciesId, bty="n", axes=T, xlab="log2(TPM)")

## intergenic
dens_intergenic <- density(log2(summed_filtered$tpm[summed_filtered$type == "intergenic"]))
dens_intergenic$y <- dens_intergenic$y * sum(summed_filtered$type == "intergenic") / length(summed_filtered$tpm)
lines(dens_intergenic, col="dodgerblue3", lwd=2)

values_G_mclust <- unique(summed_filtered$classification[summed_filtered$type == "intergenic"])
values_G_mclust <- gsub("intergenic_", "", values_G_mclust)
values_G_mclust <- as.vector(as.numeric(sort(values_G_mclust)))
summed_filtered$classification <- gsub("intergenic_","", summed_filtered$classification)

## if any point classified
for (i in 1:max(values_G_mclust)) {
  dens_intergenic_sub <- density(log2(summed_filtered$tpm[summed_filtered$type == "intergenic" & summed_filtered$classification == i]))
  ## y-axis scaling
  dens_intergenic_sub$y <- dens_intergenic_sub$y * length(summed_filtered$tpm[summed_filtered$type == "intergenic" & summed_filtered$classification == i]) / length(summed_filtered$tpm)
  lines(dens_intergenic_sub, col=cols[i], lwd=2)
  ## Print gaussian number on plot: at location of max value of gaussian
  text(dens_intergenic_sub$x[dens_intergenic_sub$y == max(dens_intergenic_sub$y)], 0.005, labels = i, col=cols[i])
}
## legend
legend("topleft", legend =  c(paste0("coding (", sum(summed_filtered$biotype %in% "protein_coding"), ")"), paste0("intergenic (", sum(summed_filtered$type == "intergenic"), ")")), lwd=2, col=c("firebrick3", "dodgerblue3"), lty=c(2,1), bty="n")
abline(v=log2(max_TPM), col="gray", lty=2, lwd=2)
dev.off()
