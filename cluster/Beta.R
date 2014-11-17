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

load("C:/Users/Ben/Dropbox/1Beta.RData")

timeF<-system.time(beta_all(comm = comm.df[1:10,],coph=cophm,traitdist=traitm)                   )

  
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