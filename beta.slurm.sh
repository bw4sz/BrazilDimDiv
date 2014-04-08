#!/bin/bash

#SBATCH -J GlobalMammals
#SBATCH -o /home1/02443/bw4sz/GlobalMammals/error.out
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -A TG-DEB130023   
#SBATCH -n 50 # Total cores

#SBATCH --mail-user=benweinstein2010@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load R_stats


echo "Begin"

ibrun Rscript < /home1/02443/bw4sz/GlobalMammals/PBDBeta.R > /home1/02443/bw4sz/GlobalMammals/PBDBeta.R

echo "End"


