#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 20
#SBATCH -p msilarge,msismall
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

regions=CBS138_s05_m03_r02_amr_genes.bed.gz
input=../bwa/Cglabrata_MEC_bwa_filtered_annotated.vcf.gz
output=Cglabrata_MEC_CGD_amr_genes.vcf

module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load htslib

bcftools view -R "${regions}" \
    "${input}" \
    |bcftools annotate \
     -a "${regions}" \
     -c CHROM,FROM,TO,GENE,- \
     -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">') \
    | bcftools annotate \
     -a "${regions}" \
     -c CHROM,FROM,TO,-,LOCUS \
     -h <(echo '##INFO=<ID=LOCUS,Number=1,Type=String,Description="Locustag">') \
    >  "${output}"

bgzip "${output}"
tabix "${output}".gz
