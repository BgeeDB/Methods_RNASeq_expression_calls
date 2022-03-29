source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")
source("./scripts/InputFiles_for_benchmark.R")

outputFolder <- "./figures/Suppl_figure_9/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

humanData <- humanData[,2:ncol(humanData)]
pdf(file = paste0(outputFolder, "/Human_lung_Density_Variance_samples.pdf"), width = 20, height = 8)
par(mfrow=c(1,2))
## plot of the density distribution
plot(density(log(humanData[,1])), main =paste0("Human lung libraries"), xlab = "log Expression", ylim=c(0,0.25), bty="n")
for(i in 1:ncol(humanData)) lines(density(log(humanData[,i])))
abline(v = c(1,4), lwd = c(2, 1), col = c("darkblue", "indianred"))
saveInfo <- c()
for (i in colnames(humanData)){
  minValue <- min(humanData[[i]][which(humanData[[i]] > 0)])
  varianceData <- var(humanData[[i]])
  keepData <- minValue < 1.e-10
  collectInfo <- c(i, minValue, varianceData, keepData)
  saveInfo <- rbind(collectInfo, saveInfo)
}
saveInfo <- as.data.frame(saveInfo)
colnames(saveInfo) <- c("Library", "Minimum_value", "Variance", "KeepLibrary")
saveInfo$Variance <- as.numeric(saveInfo$Variance)
saveInfo$Library <- as.factor(saveInfo$Library)
plot(saveInfo$Library, saveInfo$Variance, main="Human lung libraries", ylab="Variance", xlab=" ",bty="n",las = 2,cex.axis=0.7)
dev.off()
