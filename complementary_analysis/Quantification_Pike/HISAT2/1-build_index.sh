#!/bin/bash

#SBATCH --account=mrobinso_bgee
#SBATCH --partition=cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=3-00:00:00

#SBATCH --output=/Different_alignments_Pike/index_pike_hisat2.out
#SBATCH --error=/Different_alignments_Pike/index_pike_hisat2.err
#SBATCH --export=NONE
#SBATCH --job-name=index_pike_hisat2

module load gcc/9.3.0
module load python/2.7.18
module load hisat2/2.2.0

indexPath=/Different_alignments_Pike/indexes_pike
annot=/Different_alignments_Pike/InfoFolders_pike

### build index with HISAT2
hisat2-build -p 4 $annot/Esox_lucius.Eluc_v4.transcriptome.fa $indexPath/Pike_HISAT2/pike_hisat2

