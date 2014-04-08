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
require(parallel)

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
  
  #bring in traits
  traits.o <- read.table(paste(droppath,"Dimensions/Data/Brazil/BenH/All_Mammal_Data-9-6-13.txt",sep=""),header=TRUE)
  
  #Which traits do we want, just body mass and range size for now.
  traits<-traits.o[,colnames(traits.o) %in% c("logBodyMass","logAllGeogRange")]
  
  rownames(traits)<-traits.o$IUCN.binomial
  
  #just get a complete dataset
  traits<-traits[complete.cases(traits),]
  
  dim(traits)
  
  #Broadcast to all other nodes, can this be run in one command?
  #Option 1, pack them into a list and unlist them on each node
  dataExport<-list(siteXspp,tree,splist,traits)
  bcast(dataExport)
  
} else {
  x<-NULL
}

#sink output for overnight runs so we can see it tomorrow
#sink("")

###########################
###############Read in data on the node
###########################

siteXspp<-dataExport[[1]]
tree<-dataExport[[2]]
splist<-dataExport[[3]]
splist<-dataExport[[4]]

#make sure to match the tip labels and the siteXspp matrix
table(tree$tip.label %in% colnames(siteXspp))

#list loaded packages
(.packages())

###Define Source Functions
source(paste(gitpath,"BrazilSourceFunctions.R",sep=""))

##########################################################
#Read in data
##########################################################

#Do we want to subset the data for a test case? If so, add which rows below
comm<-siteXspp[1:5,]

dim(comm)


#Create index

#Get the portion of the data for that rank

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

#################################################################################
#Richness of the assemblages
richness_sites<-apply(siteXspp,1,sum)

#Composition of each site
sp.list<-apply(siteXspp,1,function(x){
  names(x[which(x==1)])
})

#Legacy correction, data.merge is data.df, sorry
data.df<-data.merge

#Write to file
write.csv(data.df,paste(droppath,"FinalData.csv",sep=""))

#Or save data
save.image(paste(droppath,"MetricsDimDiv.RData",sep=""))

#Data Generation Complete
##########################################################################################
##########################################################################################

#Depending if you are a cluster
finalize()
