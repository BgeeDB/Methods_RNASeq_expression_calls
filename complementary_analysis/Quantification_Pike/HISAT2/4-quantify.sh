#!/bin/bash

#SBATCH --account=mrobinso_bgee
#SBATCH --partition=cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=3-00:00:00

#SBATCH --output=/Different_alignments_Pike/Hisat2_Quantify.out
#SBATCH --error=//Different_alignments_Pike/Hisat2_Quantify.err
#SBATCH --export=NONE
#SBATCH --job-name=Hisat2_Quantify

module load gcc/9.3.0
module load salmon/0.12.0

### transcriptome
pathTranscriptome=/Different_alignments_Pike/InfoFolders_pike/Esox_lucius.Eluc_v4.transcriptome.fa
## path to the data
pathData=/Different_alignments_Pike/pike_folder_samples/


folders=(SRX514235 SRX514236 SRX514237 SRX514238 SRX514240 SRX514258 SRX514263 SRX514266 SRX514267 SRX514268 SRX514269 SRX514270 SRX514271 SRX667246 SRX667247 SRX667248 SRX667249 SRX667250 SRX667251 SRX667252 SRX667253 SRX667254 SRX667255 SRX667256)

for folder in "${folders[@]}"
do
cd $pathData/$folder/
echo $folder
salmon quant -t $pathTranscriptome -l A -a $pathData/$folder/align_hisat2.bam -o $pathData/$folder/Hisat2_quantify_with_salmon_quant
done

