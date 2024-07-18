#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=1:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --array=83

set -ue
set -o pipefail

line=${SLURM_ARRAY_TASK_ID}
sample_file=../Cglab_sequencing_paths.txt
tempdir=/scratch.global/scot0854/cglabrata/  # with trailing slash
species=Cglabrata
ref_fasta=/home/selmecki/shared/disaster_recovery/Reference_Genomes/Cglabrata/CBS138_s05-m03-r02/C_glabrata_CBS138_version_s05-m03-r02_chromosomes.fasta

# Read sample file line corresponding to array task ID and get variables
strain=$(awk -v val="$line" 'NR == val { print $1}' $sample_file)
read1=$(awk -v val="$line" 'NR == val { print $2}' $sample_file)
read2=$(awk -v val="$line" 'NR == val { print $3}' $sample_file)

# Load modules for trimming and aligning
module use /home/selmecki/shared/software/modulefiles.local/

module load samtools/1.10
module load bbmap

# Check for/create output directories
mkdir -p "${tempdir}"trimmed_fastq logs bbmap

# JGI BBTools data preprocessing guidelines:
## trim adapters
bbduk.sh in1="${read1}" in2="${read2}" \
out1="${tempdir}"trimmed_fastq/"${strain}"_trim_adapt1.fq \
out2="${tempdir}"trimmed_fastq/"${strain}"_trim_adapt2.fq \
ref=adapters ktrim=r k=23 mink=11 hdist=1 ftm=5 tpe tbo

## contaminant (phix) filtering per bbduk user guide
bbduk.sh in1="${tempdir}"trimmed_fastq/"${strain}"_trim_adapt1.fq \
in2="${tempdir}"trimmed_fastq/"${strain}"_trim_adapt2.fq \
out1="${tempdir}"trimmed_fastq/"${strain}"_unmatched1.fq \
out2="${tempdir}"trimmed_fastq/"${strain}"_unmatched2.fq \
outm1="${tempdir}"trimmed_fastq/"${strain}"_matched1.fq \
outm2="${tempdir}"trimmed_fastq/"${strain}"_matched2.fq \
ref=phix,artifacts k=31 hdist=1 stats=logs/"${strain}"_phistats.txt

## quality trimming (bbduk user guide recommends this as separate step from adapter trimming)
bbduk.sh in1="${tempdir}"trimmed_fastq/"${strain}"_unmatched1.fq \
in2="${tempdir}"trimmed_fastq/"${strain}"_unmatched2.fq \
out1="${tempdir}"trimmed_fastq/"${strain}"_trimmed_1P.fq \
out2="${tempdir}"trimmed_fastq/"${strain}"_trimmed_2P.fq \
qtrim=rl trimq=10

## alignment with bbmap (global aligner, may work better for longer indels)
bbmap.sh in1="${tempdir}"trimmed_fastq/"${strain}"_trimmed_1P.fq \
    in2="${tempdir}"trimmed_fastq/"${strain}"_trimmed_2P.fq \
    ref="${ref_fasta}" \
    sam=1.3 \
    rgid="${species}"_"${strain}" \
    rgsm="${strain}" \
    trd=t \
    out=bbmap/"${strain}"_bbmap.sam \
    lhist=logs/"${strain}"_lhist.txt \
    indelhist=logs/"${strain}"_indelhist.txt \
    ihist=logs/"${strain}"_ihist.txt \
    statsfile=logs/"${strain}"_stats.txt

samtools collate -@8 -O bbmap/"${strain}"_bbmap.sam \
    | samtools fixmate -@8 -m - - \
    | samtools sort -@8 -l 0 \
    | samtools markdup -@8 - bbmap/"${strain}"_bbmap_sorted_markdup.bam

samtools index bbmap/"${strain}"_bbmap_sorted_markdup.bam

rm bbmap/"${strain}"_bbmap.sam
