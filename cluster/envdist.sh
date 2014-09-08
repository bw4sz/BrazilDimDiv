#!/bin/bash

#SBATCH -J GlobalMammals
#SBATCH -o /work/02443/bw4sz/GlobalMammals/errordist2.out
#SBATCH -p development 
#SBATCH -t 02:00:00
#SBATCH -A TG-DEB130023   
#SBATCH -n 1 # Total tasks
#SBATCH -N 1 # Total tasks


#SBATCH --mail-user=benweinstein2010@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

module load Rstats

cd /work/02443/bw4sz/GlobalMammals/

echo "Begin"

#find environmental betadiversity
#Rscript env.R > dist2.out

#find geographic distance between sites
#Rscript xydist.R > xy.out

#combine datatables.
Rscript combinetables.R > combine2.out

#write giant table, it works better as a seperate session?
#ibrun Rscript writeTABLE.R > writeT.out

echo "End"