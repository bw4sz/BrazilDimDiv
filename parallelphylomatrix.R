###Converting the phylogeny into a matrix

#Ben Holt has provided alot of the context for this work.
#Everytime we compute pairwise betadiversity, we find that the function wastes alot of time

#What we want to do is run a script to get all branches by all cell matrix

##################

require(reshape)
require(picante)
require(foreach)
require(doSNOW)
require(phylobase)

##If running locally set droppath

#Set dropbox path
droppath<-"C:/Users/Ben/Dropbox/"

#Set git path
gitpath<-"C:/Users/Ben/Documents/BrazilDimDiv/"

#Set dropbox path
#droppath<-"/home1/02443/bw4sz/DimDiv/"

#Set git path
#gitpath<-"/home1/02443/bw4sz/DimDiv/"

###Define Source Functions, does this need to be run and distributed to all nodes, can they source simultaneously
source(paste(gitpath,"BrazilSourceFunctions.R",sep=""))

#Read in data on the first file

#Read in species matrix
siteXspp <- read.csv(paste(droppath,"Dimensions/Data/Brazil/BenH/ninexdat.csv",sep=""))

#Just get the species data, starts on column 33 for this example
siteXspp<-siteXspp[,33:ncol(siteXspp)]

#Remove lines with less than 2 species
richness<-apply(siteXspp,1,sum)
keep<-which(richness > 2)
siteXspp<-siteXspp[keep,]

#Get entire species list
splist<-colnames(siteXspp)

#Read in phylogeny
tree<-read.tree(paste(droppath,"Dimensions/Data/Brazil/BenH/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk",sep=""))

#remove species in the siteXspp that are not in phylogeny
siteXspp<-siteXspp[,colnames(siteXspp) %in% tree$tip.label]

#Turn phylo object into a matrix
system.time(spp_br<-phylMatrix(phyl=tree,com=siteXspp,clust=3))

dim(spp_br)


#write to file, it can't be csv to preserve headings
save(spp_br,file=paste(gitpath,"spp_br.RData",sep=""))
