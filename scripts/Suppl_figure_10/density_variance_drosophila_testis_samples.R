source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")
source("./scripts/InputFiles_for_benchmark.R")

outputFolder <- "./figures/Suppl_figure_10/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

flyData <- flyData[,2:ncol(flyData)]
pdf(paste0(outputFolder, "/Drosophila_testis_Density_Variance_samples.pdf"), width=20, height=6)
par(mfrow=c(1,2))
## plot of the density distribution
plot(density(log(flyData[,1])), main =paste0("Drosophila testis libraries"), xlab = "log Expression", ylim=c(0,0.25), bty="n")
for(i in 1:ncol(flyData)) lines(density(log(flyData[,i])))
abline(v = c(1,4), lwd = c(2, 1), col = c("darkblue", "indianred"))
saveInfo <- c()
for (i in colnames(flyData)){
  minValue <- min(flyData[[i]][which(flyData[[i]] > 0)])
  varianceData <- var(flyData[[i]])
  keepData <- minValue < 1.e-10
  collectInfo <- c(i, minValue, varianceData, keepData)
  saveInfo <- rbind(collectInfo, saveInfo)
}
saveInfo <- as.data.frame(saveInfo)
colnames(saveInfo) <- c("Library", "Minimum_value", "Variance", "KeepLibrary")
saveInfo$Variance <- as.numeric(saveInfo$Variance)
saveInfo$Library <- as.factor(saveInfo$Library)
plot(saveInfo$Library, saveInfo$Variance, main="Drosophila testis libraries", ylab="Variance", xlab=" ",bty="n",las = 2,cex.axis=0.7)
dev.off()