#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=10:00:00
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err

set -ue
set -o pipefail

# local modules

module use "$HOME"/modulefiles.local
module load flye/20230228

in_reads=/home/selmecki/shared/disaster_recovery/Sequencing_Runs/AnnaSelmecki230621/pacbio/MEC001/MEC001_pacbio.fastq.gz
out_dir=flye

flye --pacbio-hifi "${in_reads}" --out-dir "${out_dir}" --threads 4 -g 12m
