#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 45
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err

cens=$(find . -name "CBS138_cen_*.fa" | sort)

module load mafft/7.475

for c in $cens; do
    mafft --thread 1 --reorder --globalpair --maxiterate 1000 "$c" > align_"$(basename "$c")"
done
