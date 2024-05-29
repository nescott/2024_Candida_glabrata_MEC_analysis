#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msilarge,msismall
#SBATCH -o %j.out
#SBATCH -e %j.err

# Concatenate chromosome vcfs into genome-wide vcf from Freebayes variant calling
# Bcftools Filter: remove complex variants, remove fixed variants, mapping quality >40,
#  strand balance probability of alt > 0, number of reads on reverse strand > 0,
#  number of reads right/left of alt > 1
# Snpeff annotate using manually built DB
# SnpSift for high and moderate impact variants, zip and index

set -ue
set -o pipefail

# clean up intermediates
function finish {
  rm chr_files.txt
  rm samples.txt
  rm "${raw_vcf}"
  rm "${sorted_vcf}"
  rm "${bcftools_out}"*
}

trap finish EXIT

chr_dir=../chr_vcf # location of vcf files to combine
raw_vcf=Cglabrata_MEC_cgd.vcf  # name for unfiltered vcf, include .vcf
sorted_vcf=Cglabrata_MEC_sorted.vcf # name for sorted, contatenated vcf, include .vcf
bcftools_out=Cglabrata_MEC_bcf.vcf  # intermediate vcf file, include .vcf
regions=bcf_CBS138_cens.tab
cen=Cglabrata_MEC_cens.vcf

#Load modules
module load bcftools/1.10.2
module load htslib/1.9

# concatenate regions, leave uncompressed for faster piping below
find "$chr_dir" -name "*.vcf" > chr_files.txt
bcftools concat -f chr_files.txt -o "${raw_vcf}"

# sort samples alphanumerically (original freebayes sorting is random)
bcftools query -l "${raw_vcf}" | sort > samples.txt
bcftools view -S samples.txt "${raw_vcf}" > "${sorted_vcf}"

 bcftools view -i \
"INFO/MQM>=40 & INFO/SAR>=1 & INFO/SAP>0 & INFO/RPL>1 & INFO/RPR>1" \
 -o "${bcftools_out}" "${sorted_vcf}"

bgzip "${bcftools_out}"
tabix "${bcftools_out}".gz

tabix -h -R "${regions}" "${bcftools_out}".gz > "${cen}"

bgzip "${cen}"
tabix "${cen}".gz
