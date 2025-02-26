#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --output=hifiasm_%j.out
#SBATCH --error=hifiasm_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH --time=12:00:00
#SBATCH -p msismall,msilarge

module use /home/selmecki/scot0854/modulefiles.local

#Load module
module load hifiasm

hifiasm -o MEC175.asm -t 48 -l 0 --n-hap 1 --hg-size 12g \
/home/selmecki/shared/disaster_recovery/Sequencing_Runs/AnnaSelmecki230621/pacbio/MEC175/MEC175_pacbio.fastq.gz
