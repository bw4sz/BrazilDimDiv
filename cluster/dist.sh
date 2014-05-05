#!/bin/bash

#SBATCH -J GlobalMammals
#SBATCH -o /home1/02443/bw4sz/GlobalMammals/errordist.out
#SBATCH -p development 
#SBATCH -t 00:10:00
#SBATCH -A TG-DEB130023   
#SBATCH -n 5 # Total tasks

#SBATCH --mail-user=benweinstein2010@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load Rstats

cd /home1/02443/bw4sz/GlobalMammals/

echo "Begin"

ibrun Rscript commdist.R > dist.out

echo "End"


