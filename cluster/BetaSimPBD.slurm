#!/bin/bash

#SBATCH -J GlobalMammals
#SBATCH -o /home1/02443/bw4sz/GlobalMammals/error.out
#SBATCH -p normal
#SBATCH -t 2:00:00
#SBATCH -A TG-DEB130023   
#SBATCH -n 16 # Total cores

#SBATCH --mail-user=benweinstein2010@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load Rstats

cd /home1/02443/bw4sz/GlobalMammals/

echo "Begin"

Rscript BetaSimBranchesMcapply.R > BetaSimPBD.out

#just the tcellbr 
#ibrun Rscript BetaSimBranches.R > BetaSimPBD.out


#echo "Betadiversity"

#ibrun Rscript BetaSimPBD.R > BetaSimPBD.out
echo "End"


