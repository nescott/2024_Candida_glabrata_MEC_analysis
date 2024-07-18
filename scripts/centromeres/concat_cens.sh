#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=scot0854@umn.edu
#SBATCH -t 30
#SBATCH -p msismall,msilarge
#SBATCH -o %j.out
#SBATCH -e %j.err

samples=$(find . -name "*_cen.fa" | sort)
arr=("A" "B" "C" "D" "E" "F" "G" "H" "I" "J" "K" "L" "M")

for s in $samples; do
  for let in "${arr[@]}"; do
    grep -A 2 "Chr$let" "$s"  >> CBS138_cen_"$let.fa" ;
  done;
done
