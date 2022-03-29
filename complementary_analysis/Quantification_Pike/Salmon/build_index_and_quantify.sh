#!/bin/bash

#SBATCH --account=mrobinso_bgee
#SBATCH --partition=cpu

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G
#SBATCH --time=3-00:00:00

#SBATCH --output=/Different_alignments_Pike/salmon_pike.out
#SBATCH --error=/Different_alignments_Pike/salmon_pike.err
#SBATCH --export=NONE
#SBATCH --job-name=salmon_pike

module load gcc/9.3.0
module load salmon/0.12.0

### build index with salmon
pathTranscriptome=/Different_alignments_Pike/InfoFolders_pike
indexes=/Different_alignments_Pike/indexes_pike
## all samples have read > 75bp so index will be build with k=31 by default
salmon index -t $pathTranscriptome/Esox_lucius.Eluc_v4.transcriptome.fa -i $indexes/Pike_Salmon

## path to the data
pathData=/Different_alignments_Pike/pike_folder_samples/

folders=(SRX514235 SRX514236 SRX514237 SRX514238 SRX514240 SRX514258 SRX514263 SRX514266 SRX514267 SRX514268 SRX514269 SRX514270 SRX514271 SRX667246 SRX667247 SRX667248 SRX667249 SRX667250 SRX667251 SRX667252 SRX667253 SRX667254 SRX667255 SRX667256)


for folder in "${folders[@]}"
do
cd $pathData/$folder/
echo $pathData/$folder/
files=($(find . -type f -name '*.fastq.gz' -exec basename {} \;))
echo $pathData/$folder/${files[0]}
echo $pathData/$folder/${files[1]}
salmon quant -i $indexes/Pike_Salmon -l A -1 $pathData/$folder/${files[0]} -2 $pathData/$folder/${files[1]} -o $pathData/$folder/transcripts_quant_salmon
done

