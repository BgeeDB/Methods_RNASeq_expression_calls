source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

outputFolder <- "./figures/Suppl_figure_8/"
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

folderPath <- "./data/Droplet_scRNA-Seq/All_FINAL_RESULTS/"

libraryIDs <- c("SRX3815586","SRX3815587","SRX3815588","SRX3815589","SRX3815590","SRX3815591","SRX6060800","SRX6060801","SRX6060802","SRX6060803",
                "SRX6060804","SRX6060805", "SRX6060806","SRX6060807","SRX6060811","SRX6060812","SRX6060813","SRX6060814","SRX6060831","SRX6060832",
                "SRX6060833","SRX6060834","SRX6060840","SRX6060841","SRX6060842")

allInfoCPMcalls <- c()
for (library in libraryIDs) {
  
  print(library)
  callsFiles <- list.files(paste0(folderPath,library), pattern = "_genic.tsv$")
  
  speciedID <- unique(scRNAlibraryInfo$speciesId[scRNAlibraryInfo$libraryId == library])
  
  collectInfo <- c()
  for (file in callsFiles) {
    print(file)
    readcalls <- read.table(file.path(folderPath,library,file), header = TRUE, sep="\t")
    readcalls$callsCPM<- ifelse(readcalls$CPM >= 1, "present", "absent")
    
    coding <- nrow(readcalls[readcalls$biotype == "protein_coding", ])
    codingPresentUMI <- nrow(readcalls[readcalls$biotype == "protein_coding" & readcalls$callsCPM == "present", ])
    proportionPC <- (codingPresentUMI/coding)*100
    
    infoFile <- data.frame(library, speciedID, file, coding, codingPresentUMI, proportionPC)
    collectInfo <- rbind(collectInfo, infoFile)
  }
  allInfoCPMcalls <- rbind(allInfoCPMcalls, collectInfo)
}

allInfoCPMcalls$organism <- ifelse(allInfoCPMcalls$speciedID == "9606", "Homo sapiens", "Mus musculus")

pdf(file=file.path(outputFolder,"Calls_TargetBased_CPM_cutoff.pdf"),width=6, height=6)
g2 <- ggplot(allInfoCPMcalls, aes(organism, proportionPC, fill = organism)) + 
  geom_boxplot() + ylim(0,100) + xlab(" ") + ylab("% Protein coding genes present") + 
  ggtitle(expression("Calls of expressed genes with CPM">=1)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),legend.position = "none")
g2
dev.off()
