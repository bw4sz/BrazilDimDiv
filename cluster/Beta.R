#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW

##################

require(pbdMPI,quietly=TRUE,warn.conflicts=FALSE)

init()


#Require libraries
suppressMessages(require(reshape,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(picante,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(ggplot2,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(parallel,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(foreach,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(GGally,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(data.table,quietly=TRUE,warn.conflicts=FALSE))

#Everyone say hello
comm.print(comm.rank(), all.rank = TRUE)

##If running locally set droppath

droppath<-"/home1/02443/bw4sz/GlobalMammals/"

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

print(paste("Memory used:",mem()))

timeF<-system.time(beta_out<-betaPar.scatter(comm.df,Index_Space,cophm,traitm))

print(head(beta_out))

#checkpoint, write data if fails
save.image(paste("Beta/",.rank,"Beta.RData",sep=""))

#Return timing argument to console
print(timeF)

#try writing from all
comm.write(beta_out,"Output/FinalData.txt",row.names=FALSE)

#remove data file
file.remove(fil)

finalize()

###############################
#Data Generation Complete
###############################