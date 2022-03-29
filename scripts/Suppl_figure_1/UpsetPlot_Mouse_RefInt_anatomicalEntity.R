source("./scripts/checkLibraries.R")
source("./scripts/InputFiles_From_Bgee.R")

## Note: we just keep for the upset plot anatomical entities with at least 5 libraries
## cut-off applied
## Adult_mammalian_kidney --> select max 3
## Ammons_horn --> select max 3
## brain --> select max 3
## cerebellar_cortex --> select max 5
## cerebellum --> select max 1
## colon --> select max 1
## cortex_of_kidney --> select max 7
## cortical_plate --> select max 5
## duodenum --> select max 1
## embryonic_brain --> select max 3
## embryonic_facial_prominence --> select max 3
## forebrain --> select max 3
## forelimb_bud --> select max 2
## ganglionic_eminence --> (select max 8
## heart --> select max 3
## hindbrain --> select max 3
## hindlimb_bud --> select max 1
## hindlimb_stylopod_muscle --> select max 5
## intestine --> select max 4
## lens_of_camera-type_eye --> select max 3
## limb --> select max 3
## liver --> select max 3
## lung --> select max 3
## midbrain --> select max 3
## muscle_tissue --> select max 3
## neural_tube --> select max 3
## placenta --> select max 2
## primary_visual_cortex --> select max 4
## proximal_tubule --> select max 1
## skeletal_muscle_tissue --> select max 2
## spleen --> select max 3
## stomach --> select max 4
## superior_frontal_gyrus --> select max 4
## testis --> select max 3
## thymus --> select max 3
## ventricular_zone --> select max 6
## zone_of_skin --> select max 3

outputFolder <- file.path("./figures/Suppl_figure_1/")
if (!dir.exists(outputFolder)){
  dir.create(outputFolder)
} else {
  print("Already exists!")
}

pathFolder <- "./data/Mouse_Select_RefInt_based_anatomical_entity/Sum_organized_by_anatomical_entity/"

tissue_selected <- c("adult mammalian kidney", "Ammon's horn", "brain", "cerebellar cortex", "cerebellum", "colon", "cortex of kidney",
                     "cortical plate", "duodenum", "embryonic brain", "embryonic facial prominence", "forebrain", "forelimb bud", 
                     "ganglionic eminence", "heart", "hindbrain", "hindlimb bud", "hindlimb stylopod muscle", "intestine", 
                     "lens of camera-type eye", "limb", "liver", "lung", "midbrain", "muscle tissue", "neural tube", "placenta", 
                     "primary visual cortex", "proximal tubule", "skeletal muscle tissue", "spleen", "stomach", "superior frontal gyrus",
                     "testis", "thymus", "ventricular zone", "zone of skin") 
cutoff_Gaussian_selected <- c(3,3,3,5,1,1,7,5,1,3,3,3,2,8,3,3,1,5,4,3,3,3,3,3,3,3,2,4,1,2,3,4,4,3,3,6,3)
info_tissue_gaussian_cutoff <- data.frame(tissue_selected, cutoff_Gaussian_selected)


collectRefInt <- list()
for (tissue in info_tissue_gaussian_cutoff$tissue_selected) {
  
  print(tissue)
  max_cutoff <- info_tissue_gaussian_cutoff$cutoff_Gaussian_selected[info_tissue_gaussian_cutoff$tissue_selected == tissue]
  
  #tissue <- gsub("_", " ", tissue)
  readSumFile <-  read.table(paste0(pathFolder, tissue, "/sum_abundance_gene_level+fpkm+intergenic+classification_Mus musculus.tsv"), header=TRUE, sep="\t")
  select_intergenic <- dplyr::filter(readSumFile, type == "intergenic" & classification == paste0("intergenic_",max_cutoff))
  max_TPM <- max(select_intergenic$tpm)
  refIntergenic <- dplyr::filter(readSumFile, type == "intergenic", tpm <= max_TPM)
  select_intergenicRef <- list(refIntergenic[,1])
  names(select_intergenicRef)=tissue
  collectRefInt <- append(collectRefInt, select_intergenicRef)
  
}


### Ref intergenic using all samples from Bgee 15
refIn_Bgee15 <- unlist(refInt_mouse$chr_start_end)

listInput <-  list(all_libraries = refIn_Bgee15,
                   adult_mammalian_kidney = collectRefInt[[1]],
                   ammons_horn = collectRefInt[[2]],
                   brain = collectRefInt[[3]],
                   cerebellar_cortex = collectRefInt[[4]],
                   cerebellum = collectRefInt[[5]], 
                   colon = collectRefInt[[6]],
                   cortex_of_kidney = collectRefInt[[7]],
                   cortical_plate = collectRefInt[[8]], 
                   duodenum = collectRefInt[[9]], 
                   embryonic_brain = collectRefInt[[10]],
                   embryonic_facial_prominence = collectRefInt[[11]],
                   forebrain = collectRefInt[[12]],
                   forelimb_bud = collectRefInt[[13]],
                   ganglionic_eminence = collectRefInt[[14]],
                   heart = collectRefInt[[15]],
                   hindbrain = collectRefInt[[16]], 
                   hindlimb_bud = collectRefInt[[17]], 
                   hindlimb_stylopod_muscle = collectRefInt[[18]],
                   intestine = collectRefInt[[19]],
                   lens_of_camera_type_eye = collectRefInt[[20]],
                   limb = collectRefInt[[21]], 
                   liver = collectRefInt[[22]],
                   lung = collectRefInt[[23]],
                   midbrain = collectRefInt[[24]],
                   muscle_tissue = collectRefInt[[25]],
                   neural_tube = collectRefInt[[26]], 
                   placenta = collectRefInt[[27]], 
                   primary_visual_cortex = collectRefInt[[28]],
                   proximal_tubule = collectRefInt[[29]], 
                   skeletal_muscle_tissue = collectRefInt[[30]],
                   spleen = collectRefInt[[31]],
                   stomach = collectRefInt[[32]],
                   superior_frontal_gyrus = collectRefInt[[33]],
                   testis = collectRefInt[[34]],
                   thymus = collectRefInt[[35]],
                   ventricular_zone = collectRefInt[[36]],
                   zone_of_skin = collectRefInt[[37]])


pdf(file.path(outputFolder, "UpsetPlot_MusMusculus_refIntergenic.pdf"),width=25,height=15)
upset(fromList(listInput), 
      sets = c("all_libraries","adult_mammalian_kidney", "ammons_horn", "brain", "cerebellar_cortex", "cerebellum", "colon", "cortex_of_kidney",
               "cortical_plate", "duodenum", "embryonic_brain", "embryonic_facial_prominence", "forebrain", "forelimb_bud", 
               "ganglionic_eminence", "heart", "hindbrain", "hindlimb_bud", "hindlimb_stylopod_muscle", "intestine", 
               "lens_of_camera_type_eye", "limb", "liver", "lung", "midbrain", "muscle_tissue", "neural_tube", "placenta", 
               "primary_visual_cortex", "proximal_tubule", "skeletal_muscle_tissue", "spleen", "stomach", "superior_frontal_gyrus",
               "testis", "thymus", "ventricular_zone", "zone_of_skin"),
      order.by = "freq",
      query.legend = "bottom",
      sets.bar.color=c("black","black","black","black","black","black","black","black","black","black","black","black","black"
                       ,"black","black","black","black","black","black","red","black","black","black","black","black","black","black","black"
                       ,"black","black","black","black","black","black","black","black","black","black"),
      mainbar.y.label = "Intersection of reference intergenic regions", 
      sets.x.label = "Total of reference intergenic regions", 
      queries = list(
        list(
          query = intersects,
          params = list("all_libraries","adult_mammalian_kidney", "ammons_horn", "brain", "cerebellar_cortex", "cerebellum", "colon", "cortex_of_kidney",
                        "cortical_plate", "duodenum", "embryonic_brain", "embryonic_facial_prominence", "forebrain", "forelimb_bud", 
                        "ganglionic_eminence", "heart", "hindbrain", "hindlimb_bud", "hindlimb_stylopod_muscle", "intestine", 
                        "lens_of_camera_type_eye", "limb", "liver", "lung", "midbrain", "muscle_tissue", "neural_tube", "placenta", 
                        "primary_visual_cortex", "proximal_tubule", "skeletal_muscle_tissue", "spleen", "stomach", "superior_frontal_gyrus",
                        "testis", "thymus", "ventricular_zone", "zone_of_skin"), 
          color = "darkblue", 
          active = T,
          query.name = "All condictions"
        )
      )
)
dev.off()

