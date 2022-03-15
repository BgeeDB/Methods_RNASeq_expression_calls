### **Robust inference of expression state in bulk and single-cell RNA-Seq using curated intergenic regions**

Sara S. Fonseca Costa, Marta Rosikiewic, Julien Roux, Julien Wollbrett, Frederic B. Bastian, Marc Robinson-Rechavi


<img src="figures/overview.png" width="1280"/>

The paper can be found on [bioRxiv]().

#### The repository
This repository collect all the data files as well as scripts necessary to re-generate all the figures of the methods paper to call expressed genes on RNASeq data.

The repository is organized by 5 main folders:

*   [data/](data/)

Folder that contain all input data necessary to reproduce the figures of the paper

*   [figures/](figures/)

Folder that contain all the figures of the paper

*   [scripts/](scripts/)

Folder that contain all the scripts used to generate the figures and get the statistics information

*   [stats_info/](stats_info/)

Folder that contain the .tsv files with stats information 


*   [analysis_info/](analysis_info/)

Folder that contain the R version and packages version used during the analysis.


##### Raw data information

All the raw data (fastq.gz files) used on this paper are available on public databases with exception for GTEx data that can be retrieved using [BgeeDB](https://bioconductor.org/packages/release/bioc/html/BgeeDB.html) R package or through [bgee website](https://bgee.org/) but already processed.
