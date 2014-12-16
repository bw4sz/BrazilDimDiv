#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW

##################

args<-commandArgs(TRUE)

library(pbdMPI,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)

init()

#set position within the loop
position_start<-args[1]

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

paste("Running splits from",position_start," to ",position_start+comm.size())

#rank to filename
fil<-paste("Data/",.rank,"Rank.RData",sep="")

#read in portion of the data
load(fil)

#Calculate phylogenetic species lists
sp.list.phylo<-apply(comm.df[,!colnames(comm.df) %in% "id",with=F],1,function(x){ 
  names(x[which(x==1)])
})

names(sp.list.phylo)<-comm.df$id

col_keep<-colnames(comm.df)[colnames(comm.df) %in% colnames(traitm)]
comm.df.trait<-comm.df[,c(col_keep,"id"),with=F]

sp.list.trait<-apply(comm.df.trait[,!colnames(comm.df.trait) %in% "id",with=F],1,function(x){ 
  names(x[which(x==1)])
})

names(sp.list.trait)<-comm.df$id

system.time(beta_out<-betaPar.scatter(toScatterIndex = Index_Space,coph=cophm,traitdist=traitm,sp.list.phylo = sp.list.phylo,sp.list.trait = sp.list.trait))

print(head(beta_out))

#checkpoint, write data if fails
#save.image(paste("Beta/",.rank,"Beta.RData",sep=""))

#Return timing argument to console
#print(timeF)

#try writing from all
comm.write.table(beta_out,"Output/FinalDataLog1dg.txt",row.names=F,append=T)

#remove data file
file.remove(fil)

finalize()

###############################
#Data Generation Complete
###############################
