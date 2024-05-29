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
raw_vcf=Cglabrata_MEC_strains.vcf  # include .vcf
sorted_vcf=Cglabrata_sorted_MEC.vcf # include .vcf
ref=CBS138_s05m03ro2
fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/Cglabrata/CBS138_s05-m03-r02/C_glabrata_CBS138_version_s05-m03-r02_chromosomes.fasta
bcftools_out=Cglabrata_MEC_filtered.vcf  # include .vcf
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
  rm "${bcftools_out}"
 }

 trap finish EXIT

 # concatenate regions
find "$chr_dir" -name "*.vcf" | sort > chr_files.txt
bcftools concat -f chr_files.txt -o "${raw_vcf}"
bgzip "${raw_vcf}"
tabix "${raw_vcf}".gz

# sort samples alphanumerically
bcftools query -l "${raw_vcf}".gz | sort > samples.txt
bcftools view -S samples.txt "${raw_vcf}".gz \
  | bcftools view -i "INFO/MQM>39" \
  | bcftools view -i "FORMAT/DP[*] >9" \
  | bcftools view -e 'INFO/TYPE="complex"' \
  | bcftools norm -f "${fasta}" -m -indels \
  | bcftools +fill-tags - -- -t "FORMAT/VAF" \
  | bcftools view -i "FORMAT/VAF[*] > 0.9" \
  | bcftools view -i "INFO/SAR>=1 & INFO/SAP>0 & INFO/RPL>1 & INFO/RPR>1" \
  | bcftools view -c 1 \
  -o "${bcftools_out}"

# annotate using snpeff with manually built database
java -Xmx9g -jar "${snpeff}" -c "${snpeff_config}" "${snpeff_db}" \
"${bcftools_out}" >  "${annotate_vcf}"

# subset annotated to just SNPs
# and output tab-delimited file for use in R MCA script for preliminary clustering
tr "\n" "\t" < samples.txt > "${genotypes_table}"
sed -i -e '$\a' "${genotypes_table}"
sed '1s/CHROM\tPOS\t/' "${genotypes_table}"

bcftools view -e 'GT="mis"' "${annotate_vcf}" \
| bcftools view -m2 -M2 -v snps \
| bcftools query -f '%CHROM\t%POS[\t%GT]\n' >> "${genotype_table}"

bgzip "${annotate_vcf}"
tabix "${annotate_vcf}".gz
