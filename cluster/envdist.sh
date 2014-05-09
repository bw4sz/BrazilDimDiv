#!/bin/bash

#SBATCH -J GlobalMammals
#SBATCH -o /home1/02443/bw4sz/GlobalMammals/errordist.out
#SBATCH -p development 
#SBATCH -t 02:00:00
#SBATCH -A TG-DEB130023   
#SBATCH -n 20 # Total tasks
#SBATCH -N 5 # Total tasks

#SBATCH --mail-user=benweinstein2010@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load Rstats

cd /home1/02443/bw4sz/GlobalMammals/

echo "Begin"

#find environmental betadiversity
ibrun Rscript env.R > dist.out

#find geographic distance between sites
Rscript xydist.R > xy.out

#combine datatables.
Rscript combineenv.R > combine.out

echo "End"