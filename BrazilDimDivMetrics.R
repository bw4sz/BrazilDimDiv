#Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
#Manuscript submitted to the American Naturalist

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW
#Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD

#require(Rmpi)
require(doSNOW)
require(vegan)
require(reshape)
require(ape)
require(picante)
require(ggplot2)
library(snow)
library(foreach)
require(GGally)
require(stringr)
require(scales)
library(foreach)
require(raster)
require(boot)

#Are you running locally? If yes
#cluster<-makeCluster(4,"SOCK")

#Are you running on an MPI cluster? If yes
#cluster <- getMPIcluster()

# Print the hostname for each cluster member
sayhello <- function()
{
  info <- Sys.info()[c("nodename", "machine")]
  paste("Hello from", info[1], "with CPU type", info[2])
}

names <- clusterCall(cluster, sayhello)
print(unlist(names))


####################
#set paths
####################

gitpath<-""

droppath<-""

#registerDoSNOW(cluster)

###########################
###############Read in data
###########################

###Define Source Functions

#list loaded packages
(.packages())

##########################################################
#Read in data
##########################################################

#Load in the data
load(paste(droppath,"EnvCompare.RData")

#Do we want to subset the data for a test case? If so, add which rows below
comm<-siteXspp[,]

dim(comm)

ls()

######################################################
#Create a function for computing betadiversity metrics
#######################################################

#Find taxonomic betadiversity

#####################################
##Taxonomic Betadiversity
d<-as.matrix(vegdist(comm,binary=TRUE,upper=FALSE,diag=FALSE))
sorenson<-melt(d)[melt(upper.tri(d))$value,]
colnames(sorenson)<-c("To","From","Sorenson")
######################################

#####################################################
#Perform Phylogenetic and Trait Betadiversity
#####################################################

#Define Betadiversity function

beta_all<-function(comm,tree,traits){
  
  #####################################
  ##Taxonomic Betadiversity
  d<-as.matrix(vegdist(comm,binary=TRUE,upper=FALSE,diag=FALSE))
  sorenson<-melt(d)[melt(upper.tri(d))$value,]
  colnames(sorenson)<-c("To","From","Sorenson")
  ######################################
  
  #Phylosor Calculation see Bryant 2008
  phylo.matrix<-as.matrix(phylosor(comm,tree))
  diag(phylo.matrix)<-NA
  Phylosor.phylo<-melt(phylo.matrix)
  colnames(Phylosor.phylo)<-c("To","From","Phylosor.Phylo")
  Phylosor.phylo$Phylosor.Phylo<-1-Phylosor.phylo$Phylosor.Phylo
  
  #Trait frame needs to match siteXSpp table
  mon_cut<-traits[rownames(traits) %in% colnames(comm),]
  
  #There are eight species without traits, take them out for just this portion of the analysis, keep the assemblage lsit
  siteXspp_traits<-comm[,colnames(comm) %in% rownames(mon_cut)]
  
  #Zscores, standardized by sd and subtracted means
  means<-apply(mon_cut,2,mean)
  Bill<-mon_cut$Bill - means["Bill"]/sd(mon_cut$Bill)
  Mass<-mon_cut$Mass - means["Mass"]/sd(mon_cut$Mass)
  WingChord<-(mon_cut$WingChord - means["WingChord"])/sd(mon_cut$WingChord)
  z.scores<-data.frame(Bill,Mass,WingChord)
  rownames(z.scores)<-rownames(mon_cut)
  newSGdist<-dist(z.scores,method="euclidean")
  
  #If you wanted a PCA 
  #prc_traits<-prcomp(mon_cut)
  #newSGdist <- dist(prc_traits$x)
  
  source(paste(gitpath,"Sources/BenHolttraitDiversity.R")
  
  #create sp.list
  sp.list<-lapply(rownames(siteXspp_traits),function(k){
    x<-siteXspp_traits[k,]
    names(x[which(x==1)])
  })
  
  names(sp.list)<-rownames(siteXspp_traits)
  dists <- as.matrix(newSGdist)
  
  rownames(dists) <- rownames(mon_cut)
  colnames(dists) <- rownames(mon_cut)
  sgtraitMNTD <- sapply(rownames(siteXspp_traits),function(i){
    
    #Iterator count
    #print(round(which(rownames(siteXspp_traits)==i)/nrow(siteXspp_traits),3))
    
    #set iterator
    A<-i
    
    #
    out<-lapply(rownames(siteXspp_traits)[1:(which(rownames(siteXspp_traits) == i))], function(B) {MNND(A,B,sp.list=sp.list,dists=dists)})
    names(out)<-rownames(siteXspp_traits)[1:(which(rownames(siteXspp_traits) == i))]
    return(out)
  })
  
  names(sgtraitMNTD) <- rownames(siteXspp_traits)
  melt.MNTD<-melt(sgtraitMNTD)
  colnames(melt.MNTD)<-c("MNTD","To","From")
  
  #Combine with other metrics
  Allmetrics0<-merge(Phylosor.phylo,melt.MNTD,by=c("To","From"))
  Allmetrics<-merge(Allmetrics0,sorenson,by=c("To","From"))
  
  return(Allmetrics)}

system.time(beta_metrics<-beta_all(comm=comm,tree=tree,traits=mon))

#Visualizations of the beta metrics
head(beta_metrics)

#Get rid of the nestedness components
beta_metrics<-beta_metrics[,colnames(beta_metrics) %in% c("To","From","MNTD","Phylosor.Phylo","Sorenson") ]

#####################################################
#Merge Betadiversity and Environmnetal Dissimilairity
#####################################################
data.merge<-merge(compare.env,beta_metrics,by=c("To","From"))

#save for checkpointing
save.image(paste(droppath,"DimDivRevisionCluster.RData")

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
write.csv(data.df,paste(droppath,"FinalData.csv",sep="")

#Or save data
save.image(paste(droppath,"MetricsDimDiv.RData",sep="")

#Data Generation Complete
##########################################################################################
##########################################################################################

#Depending if you are a cluster
stopCluster(cluster)
#Quit
q()
#Don't save
n
#exit terminal
exit
