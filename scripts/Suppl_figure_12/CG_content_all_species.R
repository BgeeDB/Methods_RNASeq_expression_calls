source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- "./figures/Suppl_figure_12/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

pathFolder <- "./data/CG_content_Intergenic_all_species/"
species <- c(9606, 10090, 7955, 7227, 6239, 9598, 9597, 9593, 9544, 10116, 9913, 9823, 9796, 9986, 9615, 9685, 
             10141, 13616, 9258, 9031, 28377, 8364, 7237, 7240, 7740,7897, 7918, 7936, 7994, 8010, 8030, 8049, 
             8081, 8090, 8154, 8355, 9103, 9483, 9531, 9541, 9545, 9555, 9925, 9940, 9974, 10181, 30608, 32507, 52904, 60711, 69293, 105023)

collectAllInfo <- data.frame()
for (i in species) {
  
  read_RefInt <- fread(paste0(pathFolder, i, "_intergenic_CG.tsv"), header=FALSE, sep="\t")
  mean_CG_RefInt <- mean(read_RefInt$V2)
  
  read_OtherInt <- fread(paste0(pathFolder, i, "_other_intergenic_CG.tsv"), header=FALSE, sep="\t")
  mean_CG_OtherInt <- mean(read_OtherInt$V2)
  
  getInfo <- cbind(mean_CG_RefInt, mean_CG_OtherInt)
  getInfo <- as.data.frame(getInfo)
  getInfo$species <- i
  collectAllInfo <- rbind(collectAllInfo, getInfo)
}

pdf(file = file.path(outputFolder, "/CG_content_all_species.pdf"), width = 16, height = 7)
g1 <- ggplot(collectAllInfo, aes(mean_CG_RefInt, mean_CG_OtherInt, label = species)) + 
  geom_point()+ ylim(30,60)+ xlim(30,50) +
  geom_text(hjust=1.2, vjust=0)+
  ggtitle("CG content")+
  xlab("Mean of CG content in Reference Intergenic Regions")+
  ylab("Mean of CG content in Other Intergenic Regions")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                           panel.grid.minor = element_blank(), 
                                                                           axis.line = element_line(colour = "black"),legend.position = "none")
g1
dev.off()
