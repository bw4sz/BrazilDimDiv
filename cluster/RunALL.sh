#!/bin/bash

#SBATCH -J GlobalMammals
#SBATCH -o /work/02443/bw4sz/GlobalMammals/error.out
#SBATCH -p development 
#SBATCH -t 2:00:00
#SBATCH -A TG-DEB130023
#SBATCH -n 2 # Total tasks

#SBATCH --mail-user=caterina.penone@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load Rstats

cd /work/02443/bw4sz/GlobalMammals/

#echo "Create submatrix of unique comparisons"

#Rscript uniqueSxpp.R > unique.out

echo "Begin"

echo "distribute data"

#ibrun Rscript splitData.R > splitData.out

echo "Betadiversity Analysis"

#loop through and append each to the end of the document, it will take longer, but it will use less memory.

for I in 4 6
do
ibrun Rscript Beta.R $I >> Beta.out
done  


echo "End"


