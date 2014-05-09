Dimensions of Biodiversity
============

Code by Ben Weinstein
Stony Brook University
New York, USA

Aim: To compute and compare taxonomic, phylogenetic and trait betadiversity among all global grid cells for mammals.

Input Data
------------

1. A cell by species matrix is needed for species distributions

2. A phylogeny with matching taxonomy to the species distributions

3. A trait matrix with matching taxonomy to the species distributions

Code
---------

In general, these computations need to run on a supercomputing cluster.

The code below is formatted for the Stampede Supercomputing cluster at the University of Texas. All scripts and data inside the "cluster" folder were run in the stampede environment. See **RunALL.sh** for the job scheduler script.


1. uniqueSxpp.R
  * This script finds the unique comparisons among all combinations of grid cells. It creates two files, a xytable with matching unique and original grid cells, adn the lat long coord of each cell. This data will be recombined at the end
  
2. BetaSimBranches.R
  * Courtesy of Ben Holt at Imperial College London. The idea is that calculating phylogenetic similairity among taxa is slow because the phylo objects are not matrices. We first convert the entire cell by branch structure into a single matrix, and then compute the betadiversity with the next script. This saves alot of time.
  
3. scatterBeta
  * The main workhourse function, which calls a series of source scripts from Input/BrazilSourceFunctions.R. The unique cell table is broken into pieces, based on the number of cores available, and each core computes the taxonomic, phylogenetic and trait betadiversity combinations on its own piece. The results are recombined with the original xytable and  written with their unique code and original cell number to finaldata.csv

4. envdist.R

  * Perform a pca the environmental space. The distance function is much faster, so we can just parallelize across a few cores. FinalData is recombined by environmental betadiversity based on the xycoordinates. 

**Export the data from the cluster to your local machine**

Run locally using the amazing data.table package we can compute a number of statistics. See stack overflow and the data.table FAQ, which i just learned about last week for examples. *Its worth learning. *
4.Rawresults.md file for an example of the power of data.table and some thought about potential analysis. 

