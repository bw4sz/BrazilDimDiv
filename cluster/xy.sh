#!/bin/bash

#SBATCH -J GlobalMammals
#SBATCH -o /home1/02443/bw4sz/GlobalMammals/error.out
#SBATCH -p development 
#SBATCH -t 1:00:00
#SBATCH -A TG-DEB130023
#SBATCH -n 1 # Total tasks

#SBATCH --mail-user=benweinstein2010@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load Rstats

cd /home1/02443/bw4sz/GlobalMammals/

echo "Begin"

#Rscript xydist.R > xy.out

Rscript combineenv.R > combine.out

echo "end"