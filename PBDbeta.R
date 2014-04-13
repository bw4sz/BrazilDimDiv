#Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
#Manuscript submitted to the American Naturalist

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW
#Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD


#############
#Set testing flag, 
#this subsitutes the real data for a tiny dummy example for debugging

testing<-TRUE

##################

require(pbdMPI)

init()

require(vegan)
require(reshape)
require(picante)
require(ggplot2)
require(stringr)
require(scales)
require(raster)
require(parallel)
require(foreach)

#Everyone say hello
#comm.print(comm.rank(), all.rank = TRUE)

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

#if (comm.rank()==0){ # only read on process 0
  ################
  if(testing=FALSE){
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
  
  #bring in traits
  traits.o <- read.table(paste(droppath,"Dimensions/Data/Brazil/BenH/All_Mammal_Data-9-6-13.txt",sep=""),header=TRUE)
  
  #Which traits do we want, just body mass and range size for now.
  traits<-traits.o[,colnames(traits.o) %in% c("logBodyMass","logAllGeogRange")]
  
  rownames(traits)<-traits.o$IUCN.binomial
  
  #just get a complete dataset
  traits<-traits[complete.cases(traits),]
  }
  
  #Testing data
  if(testing==TRUE)    {
    #Dummy phylogeny, siteXspp and traits from phylocom picante package
    data(phylocom)
    tree<-phylocom$phylo
    comm<-phylocom$sample
    traits<-phylocom$traits
    print("Testing Data Loaded")
  }
 
  #################Data prep on first rank, we want to place the most time consuming steps that are not parallelizable here?
  #Turn phylo object into a matrix
  system.time(branch.out<-branching(tree))
  
  system.time(tcellbr<-psimbranches(tree,comm,branch.out))
  
  ####Data Load Complete####
  #compute cell by branch relationship for all cells for ben holt's function
  
  #Broadcast to all other nodes, can this be run in one command?
  #Option 1, pack them into a list and unlist them on each node
  dataExport<-list(siteXspp,tree,splist,traits,tcellbr)
  
} else {
  dataExport<-NULL
}

#broadcast to all files
bcast(dataExport)

#sink output for overnight runs so we can see it tomorrow
#sink("")

########################################
###############Read in data on the node
########################################

#Seperate and name the objects that were broadcast from the rank 0 node.
siteXspp<-dataExport[[1]]
tree<-dataExport[[2]]
splist<-dataExport[[3]]
traits<-dataExport[[4]]
tcellbr<-dataExport[[5]]
#list loaded packages
(.packages())


#Do we want to subset the data for a test case? 
#If so, add which rows below
comm<-siteXspp[1:20,]

#####################################################
##############Compute Betadiversity##################
#####################################################

#betaPar takes in three arguments, the siteXspp matrix, the rank number of the node, and the number of chunks
#Until we can get on the cluster the 2nd argument will be 1 for testing
#rank<-comm.rank
#system.time(beta_out<-betaPar(siteXspp,rank,2))


timeF<-system.time(beta_out<-betaPar(rbind(comm,comm),1,2,beta.sim=TRUE,phylosor.c=FALSE))
timeF
timeF<-system.time(beta_out<-betaPar(rbind(comm,comm),1,2,beta.sim=FALSE,phylosor.c=TRUE))
timeF

#benchmark the two approaches across 10 runs.
require(rbenchmark)
a<-benchmark(replications=2,
          betaPar(rbind(comm,comm),1,2,beta.sim=TRUE,phylosor.c=FALSE),
          betaPar(rbind(comm,comm),1,2,beta.sim=FALSE,phylosor.c=TRUE))

#Correlation among dimensions
ggpairs(beta_out[,-c(1,2)],cor=TRUE)

#profile
Rprof(tmp <- tempfile())
beta_out<-betaPar(comm,1,2,beta.sim=TRUE,phylosor.c=TRUE)
Rprof(NULL)
summaryRprof(tmp)
unlink(tmp)


#Return timing argument to console
print(timeF)

#Gather all to comm Rank 0 and write to file? cast into giant dataframe?

if (comm.rank()==0){ # only read on process 0
  
  #######Gather from all runs#########
  beta_gather<-gather()
  all_out<-rbind.fill.matrix(beta_gather)
  write.csv(data.df,paste(droppath,"FinalData.csv",sep=""))
} else {
  x<-NULL
}

#Data Generation Complete
##########################################################################################
##########################################################################################

finalize()
