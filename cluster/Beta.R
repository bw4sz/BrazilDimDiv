#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW

##################

library(pbdMPI,quietly=TRUE,warn.conflicts=FALSE)

init()


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
.rank<-comm.rank()+1

#rank to filename
fil<-paste("Data/",.rank,"Rank.RData",sep="")

#read in portion of the data
load(fil)

#Calculate phylogenetic species lists
sp.list.phylo<-apply(comm.df,1,function(x){ 
  names(x[which(x==1)])
})

comm.df.trait<-comm.df[,colnames(comm.df) %in% rownames(traitdistance)]

sp.list.trait<-apply(comm.df.trait,1,function(x){ 
  names(x[which(x==1)])
})

system.time(beta_out<-betaPar.scatter(toScatterMatrix = comm.df,Index_Space[,1:20],coph=cophm,traitdist=traitm,sp.list.phylo = sp.list.phylo,sp.list.trait = sp.list.trait))

print(head(beta_out))

#checkpoint, write data if fails
save.image(paste("Beta/",.rank,"Beta.RData",sep=""))

#Return timing argument to console
print(timeF)

#try writing from all
comm.write.table(beta_out,"Output/FinalData.txt",row.names=F)

#remove data file
file.remove(fil)

finalize()

###############################
#Data Generation Complete
###############################