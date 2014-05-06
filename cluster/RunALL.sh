#!/bin/bash

#SBATCH -J GlobalMammals
#SBATCH -o /home1/02443/bw4sz/GlobalMammals/error.out
#SBATCH -p development 
#SBATCH -t 2:00:00
#SBATCH -A TG-DEB130023
#SBATCH -n 100 # Total tasks

#SBATCH --mail-user=benweinstein2010@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load Rstats

cd /home1/02443/bw4sz/GlobalMammals/

echo "Begin"

echo "Create submatrix of unqiue comparisons"

#Rscript uniqueSxpp.R > unique.out

echo "Find tcellbr from branches"

#ibrun  RMPISNOW < BetaSimBranches.R > BetaSimPBD.out

echo "distribute data"
ibrun Rscript splitData.R > splitData.out

echo "Betadiversity Analysis"

ibrun Rscript Beta.R > Beta.out

echo "Aggregate"

ibrun Rscript AggregateMean.R > Mean.out

echo "End"


