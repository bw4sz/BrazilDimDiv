#Find unique richnesses in the matrix

#Require libraries
suppressMessages(require(data.table,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(picante,quietly=TRUE,warn.conflicts=FALSE))

##set path
#droppath<-"/work/02443/bw4sz/GlobalMammals/"
droppath<-"C:/Users/Ben/Documents/BrazilDimDiv/cluster/"
setwd(droppath)

########################################
###############Read in data on the node 0
########################################

##Read siteXspp table
siteXspp <- fread("Input/siteXspp1dgr.csv")

#make first identification rowname
siteXspp[,V1:=1:nrow(siteXspp)]

#match with phylogeny
tree<-read.tree("Input/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")
siteXspp<-siteXspp[,colnames(siteXspp) %in% c(tree$tip.label,"V1"),with=F]

#calculate richness
siteXspp$rich<-rowSums(siteXspp[,!colnames(siteXspp) %in% "V1",with=F])

#setkey to richness
setkey(siteXspp,rich)

#Remove lines with less than 2 species
siteXspp<-siteXspp[rich>1]

#Keep only original row (V1) and richness
siteXrich<-siteXspp[,colnames(siteXspp) %in% c("V1","rich"),with=F]

#write siteXrich table
write.csv(siteXrich,"Output/siteXrich.csv")

#remove duplicates
dt.unique<-unique(siteXrich)

#Get entire species list
splist<-colnames(siteXspp[, c("x","y","rich","V1"):=NULL,])

rm(siteXspp)

##We need a null model from every richness type to every other richness type. Including to itself.
ident<-expand.grid(dt.unique$rich,dt.unique$rich)
dim(ident)
#should be roughly
200* 200

#looks good.

#Okay so now build a function that randomly draws that many species for each community and computes the metric
x<-1

sampleS<-function(x,ident){
  row<-ident[x,]
  #Get species for assemblage A
  A<-sample(splist,size=row[[1]],replace=F)
  # Get species for Assemblage B
  B<-sample(splist,size=row[[2]],replace=F)
  
  #Format the data frame
}

library(pbdMPI,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)

init()

#set position within the loop


#library libraries
suppressMessages(library(reshape2,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(picante,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(ggplot2,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(parallel,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(foreach,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(GGally,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(data.table,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(betapart,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(plyr,quietly=TRUE,warn.conflicts=FALSE))

#Everyone say hello
comm.print(comm.rank(), all.rank = TRUE)

##If running locally set droppath

droppath<-"/work/02443/bw4sz/GlobalMammals/"

setwd(droppath)

###Define Source Functions
suppressMessages(source("Input/BrazilSourceFunctions.R"))

comm.print(paste("Total nodes is:", comm.size()))

############ Read in Data

#set rank
.rank<-comm.rank()+(1+position_start)
