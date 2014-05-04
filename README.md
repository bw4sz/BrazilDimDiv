BrazilDimDiv
============

Comparing taxonomic trait and phylogenetic dissimilarity among mammals and environment
Dimensions of Biodiversity - A colloboration between the Costa and Graham labs.


Code by Ben Weinstein
Stony Brook University
New York, USA

Input Data
------------

A cell by species matrix is needed for species distributions

A phylogeny with matching taxonomy to the species distributions

A trait matrix with matching taxonomy to the species distributions

Code
---------

* uniquesSxpp.R 
** Takes in the cell by species matrix and finds all unique rows using the data.table fast reading. We then create a key and a xytable for later matching
