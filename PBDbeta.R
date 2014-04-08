#Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
#Manuscript submitted to the American Naturalist

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW
#Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD

require(vegan)
require(reshape)
require(picante)
require(ggplot2)
require(stringr)
require(scales)
require(raster)
require(boot)
require(pbdMPI)

init()

#Everyone say hello
comm.print(.comm.rank, all.rank = TRUE)

#Set dropbox path
droppath<-"C:/Users/Ben/Dropbox/"

#Set git path
gitpath<-"C:/Users/Ben/Documents/BrazilDimDiv/"

#Read in data on the first file

if (comm.rank()==0){ # only read on process 0
  ################
  
  #Read in species matrix
  siteXspp <- read.csv(paste(droppath,"Dimensions/Data/Brazil/Mammal Distribution/BrazilSiteXSpp_longlat.csv",sep=""))
  
  #Get entire species list
  splist<-colnames(siteXspp)
  
  #Read in phylogeny
  tree<-read.tree(paste(droppath,"Dimensions/Data/Brazil/New_Mammal_Timetree_IUCN_2Feb2013.nwk",sep=""))
  
  #bring in traits
  mon <- read.csv("")
  
  #Broadcast to all other nodes, can this be run in one command?
  bcast(siteXspp)
  bcast(tree)
  bcast(splist)
} else {
  x<-NULL
}



#sink output for overnight runs so we can see it tomorrow
#sink("")

###########################
###############Read in data
###########################

#make sure to match the tip labels and the siteXspp matrix
table(tree$tiplabels %in% colnames(siteXspp))

#list loaded packages
(.packages())

#load(paste(droppath,"EnvCompare.RData",sep=""))

###Define Source Functions
source(paste(gitpath,"BrazilSourceFunctions.R",sep=""))

##########################################################
#Read in data
##########################################################

#Do we want to subset the data for a test case? If so, add which rows below
comm<-siteXspp[1:5,]

dim(comm)

ls()


#Create index


#Get the portion of the data for that rank


######################################################
#Create a function for computing betadiversity metrics
######################################################

#####################################################
#Compute Betadiversity
#####################################################

system.time(beta_metrics<-beta_all(comm=comm,tree=tree,traits=mon))

#Visualizations of the beta metrics
head(beta_metrics)

#Get rid of the nestedness components
beta_metrics<-beta_metrics[,colnames(beta_metrics) %in% c("To","From","MNTD","Phylosor.Phylo","Sorenson") ]

#####################################################
#Merge Betadiversity and Environmnetal Dissimilairity
#####################################################
data.merge<-merge(compare.env,beta_metrics,by=c("To","From"))

###
#End Part 1
###

#Begin Part 2
load("/home1/02443/bw4sz/DimDiv/DimDivRevisionCluster.RData")
#################################################################################
#Richness of the assemblages
richness_sites<-apply(siteXspp,1,sum)

#Composition of each site
sp.list<-apply(siteXspp,1,function(x){
  names(x[which(x==1)])
})

#Legacy correction, data.merge is data.df, sorry
data.df<-data.merge

setwd("/home1/02443/bw4sz/")

#Write to file
write.csv(data.df,paste(droppath,"FinalData.csv",sep=""))

#Or save data
save.image(paste(droppath,"MetricsDimDiv.RData",sep=""))

#Data Generation Complete
##########################################################################################
##########################################################################################

#Depending if you are a cluster
finalize()
