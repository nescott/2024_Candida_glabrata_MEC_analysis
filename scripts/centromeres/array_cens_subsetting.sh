#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 30
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=1-98

set -ue
set -o pipefail

module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load samtools/1.14

ref_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/Cglabrata/CBS138_s05-m03-r02/C_glabrata_CBS138_version_s05-m03-r02_chromosomes.fasta.gz
region_file=CBS138_s05_m03_r02_cens.txt
sample_file=../../cglab_todo.txt
line=${SLURM_ARRAY_TASK_ID}
vcf=Cglabrata_MEC_cens.vcf.gz

#Read sample file line corresponding to array task ID and get variables
sample=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)

samtools faidx "${ref_fasta}" -r "${region_file}" \
    | bcftools consensus -s "${sample}" -I "${vcf}" > "${sample}"_cen.fa

sed -i "s/^>/>$sample\_/g" "${sample}"_cen.fa
