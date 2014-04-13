###Converting the phylogeny into a matrix

#Ben Holt has provided alot of the context for this work.
#Everytime we compute pairwise betadiversity, we find that the function wastes alot of time

#What we want to do is run a script to get all branches by all cell matrix

##################
require(pbdMPI)

init()

require(reshape)
require(picante)
require(foreach)
require(phylobase)


#################Data prep on first rank, we want to place the most time consuming steps that are not parallelizable here?
#Turn phylo object into a matrix
system.time(branch.out<-branching(tree))

#parallel R using foreach and snow

system.time(tcellbr<-psimbranches.pbd(tree,siteXspp,branch.out,))

#write to file
write.csv(tcellbr,paste(gitpath,"tcellbr.csv",sep=""))

q()
n

write