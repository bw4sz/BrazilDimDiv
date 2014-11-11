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

4. A cophenetic matrix of the phylogeny.

Code
---------

In general, these computations need to run on a supercomputing cluster.

The code below is formatted for the Stampede Supercomputing cluster at the University of Texas. All scripts and data inside the "cluster" folder were run in the stampede environment. See **RunALL.sh** for the job scheduler script.


1. uniqueSxpp.R
  * This script finds the unique comparisons among all combinations of grid cells. It creates two files, a xytable with matching unique and original grid cells, adn the lat long coord of each cell. This data will be recombined at the end

2. splitData.R

To acheive massive parallelization, we can split the ~ 78 million comparisons into chunks to be run seperately and recombined at the end. This is done by breaking the community matrix into pieces and passing the correct species traits and relatedness matrices to each node. 

3. Beta.R
  * The main workhourse function, which calls a series of source scripts from Input/BrazilSourceFunctions.R. The unique cell table is broken into pieces, based on the number of cores available, and each core computes the taxonomic, phylogenetic and trait betadiversity combinations on its own piece. The results are recombined with the original xytable and  written with their unique code and original cell number to finaldata.csv


A second slurm script is used to combine the betadiversity calculations with the distance and env calculations
see **envdist.sh**

This calls env.R, xydist.R and combinetables.R


4. env.R
 * Computes pairwise env betadiversity on the 19 biocliom variables for all combinations of cells

5. xydist.R
 * Computes distance on the earth's surface in the molleweide projection for all combinations of cells

6. Combinetables.R
 *Combines the betadiversity, distance and env tables and saves the output as a RData object. 


**Export the data from the cluster to your local machine**

Run locally using the amazing data.table package we can compute a number of statistics. See stack overflow and the data.table FAQ, which i just learned about last week for examples. *Its worth learning. *
4.Rawresults.md file for an example of the power of data.table and some thought about potential analysis. 


License: This code is available under a BSD 2-Clause License.

Copyright (c) 2013. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Contact information Ben Weinstein's email: benweinstein2010@gmail.com
