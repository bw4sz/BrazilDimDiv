#!/bin/bash

#SBATCH -J GlobalMammals
#SBATCH -o /home1/02443/bw4sz/GlobalMammals/error.out
#SBATCH -p normal
#SBATCH -t 10:00:00
#SBATCH -A TG-DEB130023   
#SBATCH -n 500 # Total cores

#SBATCH --mail-user=benweinstein2010@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load Rstats

cd /home1/02443/bw4sz/GlobalMammals/

echo "Begin"
echo "Find tcellbr from branches"

#ibrun  RMPISNOW < BetaSimBranches.R > BetaSimPBD.out

#echo "Betadiversity Analysis"

ibrun Rscript RunALL.R > PBDall.out

#Compute stats
ibrun Rscript AggregateMean.R > Mean.out

echo "End"


