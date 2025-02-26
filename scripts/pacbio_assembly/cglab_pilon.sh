#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24gb
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msismall,msilarge

#Load module
module load pilon/1.22

java -Xmx24G -jar /panfs/roc/msisoft/pilon/1.22/bin/pilon-1.22.jar \
--genome mec001_flye_pilon1.fasta \
--frags ../bam/MEC001_trimmed_bwa_sorted_markdup.bam \
--threads 16 \
--output mec001_flye_pilon2 \
--fix all
