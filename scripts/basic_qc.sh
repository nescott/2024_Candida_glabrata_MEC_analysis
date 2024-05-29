#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=4gb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %x_%u_%j.out
#SBATCH -e %x_%u_%j.err

set -ue
set -o pipefail
unset DISPLAY

file_name="$(date +%F)_Cglabrata_qc"
temp_dir=/scratch.global/scot0854

# Use local modules
module use /home/selmecki/scot0854/modulefiles.local

# Load modules
#module load qualimap/20221111
module load multiqc/20221111

#find bam '(' -name "MEC084*.bam" -o -name "MEC086*.bam" -o -name "MEC096*.bam" ')' \
#-exec qualimap bamqc -bam {} -outdir $temp_dir/{} \
#--java-mem-size=4G \;

#find bam '(' -name "MEC111*.bam" -o -name "MEC136*.bam" -o -name "MEC138*.bam" ')' \

#-exec qualimap bamqc -bam {} -outdir $temp_dir/{} \
#--java-mem-size=4G \;

#find bam '(' -name "MEC158*.bam" -o -name "MEC014*.bam" ')' \
#-exec qualimap bamqc -bam {} -outdir $temp_dir/{} \
#--java-mem-size=4G \;

# multiqc

multiqc "$temp_dir"/bam "$temp_dir"/qc_and_classifier logs/ -n "$file_name"
