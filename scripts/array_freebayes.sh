#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=22:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --array=13

# call variants for all samples in a population using freebayes, subsetting by region
set -ue
set -o pipefail

species=Cglabrata
ref=CBS138_s05m03r02
ref_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/Cglabrata/CBS138_s05-m03-r02/C_glabrata_CBS138_version_s05-m03-r02_chromosomes.fasta
line=${SLURM_ARRAY_TASK_ID}
bam_list=bam.files
region_list=Cglabrata_chroms.txt

chr=$(awk -v val="$line" 'NR == val { print $0}' $region_list)

mkdir -p chr_vcf

#Load modules
module load samtools/1.10
module load freebayes/20180409

freebayes -f "${ref_fasta}" \
  -C 10 \
  -F 0.9 \
  -p 1 \
  -r "${chr}" \
  -L "${bam_list}" \
  -v chr_vcf/"${species}_${ref}_${chr:0:11}.vcf"
