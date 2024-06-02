#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=2:00:00
#SBATCH -p amdsmall,amdlarge,amd512
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

export BCFTOOLS_PLUGINS=/home/selmecki/shared/software/software.install/bcftools/1.17/plugins

chr_dir=chr_vcf # location of vcf files to combine
raw_vcf=Cglabrata_MEC_strains.vcf.gz
gene_file=Cglabrata_CBS138_s05m03r02_genes.bed.gz
gene_vcf=Cglabrata_MEC_info_strains.vcf.gz
sorted_vcf=Cglabrata_sorted_MEC.vcf.gz
ref=CBS138_s05m03ro2
fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/Cglabrata/CBS138_s05-m03-r02/C_glabrata_CBS138_version_s05-m03-r02_chromosomes.fasta
bcftools_out=Cglabrata_MEC_filtered.vcf.gz
snpeff=/home/selmecki/shared/software/snpEff/snpEff.jar
snpeff_config=/home/selmecki/shared/software/snpEff/snpEff.config
snpeff_db=CBS138_s05m03r02
annotate_vcf=Cglabrata_MEC_bwa_filtered_annotated.vcf # include .vcf
genotype_table=Cglabrata_MEC_bwa_SNPs.txt  # tab delimited table for use with R scripts (MCA using FactoMineR)

# Load modules
module use /home/selmecki/shared/software/modulefiles.local

module load bcftools/1.17
module load htslib/1.9

# remove intermediate files
function finish {
  rm samples.txt
  rm chr_files.txt
  rm "${raw_vcf}"*
  rm "${gene_vcf}"*
  rm "${sorted_vcf}"*
  rm "${bcftools_out}"*
 }

 trap finish EXIT

 # concatenate regions
find "$chr_dir" -name "*.vcf" | sort > chr_files.txt
bcftools concat -f chr_files.txt -Oz -o "${raw_vcf}"
bcftools index "${raw_vcf}"

# Add gene name annotation as info column.
bcftools annotate \
    -a "${gene_file}" \
    -c CHROM,FROM,TO,GENE \
    -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">')\
    -o "${gene_vcf}" \
    "${raw_vcf}"

bcftools index "${gene_vcf}"

# sort samples alphanumerically
bcftools query -l "${gene_vcf}" | sort > samples.txt
bcftools view -S samples.txt -Ou "${gene_vcf}" \
  | bcftools view -i "INFO/MQM>39" -Ou \
  | bcftools view -i "FORMAT/DP[*] >9" -Ou \
  | bcftools view -e 'INFO/TYPE="complex"' -Ou \
  | bcftools norm -f "${fasta}" -m -indels -Ou \
  | bcftools +fill-tags - -- -t "FORMAT/VAF" \
  | bcftools view -i "FORMAT/VAF[*] > 0.9" -Ou \
  | bcftools view -i "INFO/SAR>=1 & INFO/SAP>0 & INFO/RPL>1 & INFO/RPR>1" -Ou \
  | bcftools view -c 1 \
  -Oz -o "${bcftools_out}"

bcftools index "${bcftools_out}"

# annotate using snpeff with manually built database
java -Xmx9g -jar "${snpeff}" -c "${snpeff_config}" "${snpeff_db}" \
"${bcftools_out}" >  "${annotate_vcf}"

# subset annotated to just SNPs
# and output tab-delimited file for use in R MCA script for preliminary clustering
module unload bcftools
module load bcftools/1.9

bcftools view -m2 -M2 -v snps "${annotate_vcf}" \
| bcftools view -e 'GT="mis"' \
| bcftools query -H -f '%CHROM\t%POS[\t%GT]\n' > "${genotype_table}"

sed -i '1s/\[[0-9]\+\]//g' "${genotype_table}"
sed -i '1s/\:GT//g' "${genotype_table}"
sed -i '1s/^\# //' "${genotype_table}"

bgzip "${annotate_vcf}"
tabix "${annotate_vcf}".gz
