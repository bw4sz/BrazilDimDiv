#Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
#Manuscript submitted to the American Naturalist

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW
#Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD


#############
#Set testing flag, 
#this subsitutes the real data for a tiny dummy example for debugging

testing<-FALSE

##################

require(pbdMPI,quietly=TRUE,warn.conflicts=FALSE)

init()

require(reshape,quietly=TRUE,warn.conflicts=FALSE)
require(picante,quietly=TRUE,warn.conflicts=FALSE)
require(ggplot2,quietly=TRUE,warn.conflicts=FALSE)
require(parallel,quietly=TRUE,warn.conflicts=FALSE)
require(foreach,quietly=TRUE,warn.conflicts=FALSE)
require(GGally,quietly=TRUE,warn.conflicts=FALSE)

#Everyone say hello
comm.print(comm.rank(), all.rank = TRUE)

##If running locally set droppath

droppath<-"/home1/02443/bw4sz/GlobalMammals/"

setwd(droppath)

###Define Source Functions, does this need to be run and distributed to all nodes, can they source simultaneously

source("Input/BrazilSourceFunctions.R")

#Read in data on the first file

if (comm.rank()==0){ # only read on process 0
  ################
  if(testing==FALSE){
  #Read in species matrix
  siteXspp <- read.csv("Input/siteXspp1dgr.csv")
  
  #xytable
xytable<-data.frame(rownames(siteXspp),siteXspp$x,siteXspp$y)

  #Just get the species data, starts on column 33 for this example
  siteXspp<-siteXspp[,!colnames(siteXspp) %in% c("X","x","y","rich")]

  #Remove lines with less than 2 species
  richness<-apply(siteXspp,1,sum)
  keep<-which(richness > 1)
  siteXspp<-siteXspp[keep,]
  
  #Get entire species list
  splist<-colnames(siteXspp)
  
  #Read in phylogeny
  tree<-read.tree("Input/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")
  
  #Read in cell by branch table made from source
  load(file="Input/tcellbr.RData")
  
  #remove species in the siteXspp that are not in phylogeny
  siteXspp<-siteXspp[,colnames(siteXspp) %in% tree$tip.label]
  
  #bring in traits
  traits <- read.table("Input/imputedmammals28apr14.txt",header=TRUE,row.names=1)
  
head(traits)
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
   
  ####Data Load Complete####
  #compute cell by branch relationship for all cells for ben holt's function
  
  #Broadcast to all other nodes, can this be run in one command?
  #Option 1, pack them into a list and unlist them on each node
  dataExport<-list(siteXspp,tree,splist,traits,tcellbr,xytable)
  
} else {
  dataExport<-NULL
}

#broadcast to all files
dataExport<-bcast(dataExport)


########################################
###############Read in data on the node
########################################
# 
# #Seperate and name the objects that were broadcast from the rank 0 node.
siteXspp<-dataExport[[1]]
tree<-dataExport[[2]]
splist<-dataExport[[3]]
traits<-dataExport[[4]]
tcellbr<-dataExport[[5]]
xytable<-dataExport[[6]]

#list loaded packages
# (.packages())
 

#Do we want to subset the data for a test case? 
#If so, add which rows below
comm<-siteXspp[1:1000,]

#####################################################
##############Compute Betadiversity##################
#####################################################

comm.print(paste("Total nodes is:", comm.size()))

.rank<-comm.rank() + 1
timeF<-system.time(beta_out<-betaPar(comm,rankNumber=.rank,chunks=comm.size(),phylosor.c=FALSE,beta.sim=TRUE))

#combine rownames with xycoordinates


#Test alterantive methods
#benchmark the two approaches across 10 runs.
#require(rbenchmark)
#a<-benchmark(replications=1,
          #betaPar(comm,1,2,beta.sim=TRUE,phylosor.c=FALSE),
          #betaPar(comm,1,2,beta.sim=FALSE,phylosor.c=TRUE))

#Correlation among dimensions
#ggpairs(beta_out[,-c(1,2)],cor=TRUE)

#Return timing argument to console
print(timeF)

#try writing from all
comm.write.table(beta_out,"beta_out.csv")

finalize()

#Data Generation Complete
##########################################################################################
##########################################################################################

