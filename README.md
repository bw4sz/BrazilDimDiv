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
  * This script finds the unique comparisons among all combinations of grid cells. It creates two files, a xytable with matching unique and original grid cells, adn the lat long coord of each cell. This data will be recombined at the end (see combinetable.R)
  
2. BetaSimBranches.R
  * Courtesy of Ben Holt at Imperial College London. The idea is that calculating phylogenetic similairity among taxa is slow because the phylo objects are not matrices. We first convert the entire cell by branch structure into a single matrix, and then compute the betadiversity with the next script. This saves alot of time.

3. Beta.R
  * The main workhourse function, which calls a series of source scripts from Input/BrazilSourceFunctions.R. The unique cell table is broken into pieces, based on the number of cores available, and each core computes the taxonomic, phylogenetic and trait betadiversity combinations on its own piece. The results are written with their unique code to beta_out.csv
  
4. Combinetable.R
  * At this point we have all betadiveristy combinations, but we need to match them back to the original cells. We use the xytable created in the first script to merge the data.table structures.

5. AggregateMean.R
  * We compute basic mean, var, and quantile statistic for each cell and return a file finaldata.csv to send to coauthors.
