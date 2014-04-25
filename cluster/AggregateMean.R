#Compute average betadiversity for each cell, both in terms of mean value and quantiles for each pair. 

require(pbdMPI)

init()
require(parallel)
require(reshape2)

#read in data

droppath<-"/home1/02443/bw4sz/GlobalMammals/"

#set rank number and size.
.rank<-comm.rank()
.size<-comm.size()

#read in data on rank 0

if(.rank==0){

#read in beta file from the cluster
dat5<-read.table(paste(droppath,"beta_out.csv",sep="/"),nrows=5)
classes <- sapply(dat5, class)

dat<-read.table(paste(droppath,"beta_out.csv",sep="/"),colClasses=classes,row.names=NULL)

#create a list of all unique cells from both columns
uniA<-unique(dat$To)

uniB<-unique(dat$From)

#find unique in both
cells<-unique(uniA,uniB)

tobcast<-list(dat,cells)
} else{
tobcast<-NULL
}

#broadcast to all nodes
tobcast<-bcast(tobcast)

#unpack the data
dat<-tobcast[[1]]
cells<-tobcast[[2]]

#remove from memory
rm(tobcast)

############Proceed With Calculations#############

#chunk list of cells by index

indexS<-splitIndices(length(cells),.size)

#have the cluster pull in space based on ranknumber
comm.index<-indexS[[.rank]]

#loop through each in the index

meanC<-sapply(comm.index,function(x){
  
  #find cell value
  target<-cells[x]
  
  #find mean
  cellX<-dat[dat$To %in% target | dat$From %in% target,]
  meanCell<-sapply(cellX[,c("Sorenson","BetaSim","MNTD")],mean,na.rm=TRUE)
  
})

colnames(meanC)<-cells[comm.index]

comm.write.table(meanC,"meanCelltable.csv",row.names=FALSE)

##############################################
#Compute quantile for each combination of cells
##############################################

#for each row in data create a list position

indexAll<-splitIndices(nrow(dat),.size)

#have the cluster pull in based on rank number
comm.indexAll<-indexAll[[.rank]]

ecTrait<-ecdf(dat$MNTD)
ecPhylo<-ecdf(dat$BetaSim)
ecTax<-ecdf(dat$Sorenson)

quantframe<-lapply(comm.indexAll,function(x){
  target<-dat[x,]
  traitQ<-ecTrait(target$MNTD)
  phyloQ<-ecTrait(target$BetaSim)
  taxQ<-ecTrait(target$Sorenson)
  return(data.frame(To=target$To,From=target$From,traitQ,phyloQ,taxQ))
})

beta_quantiles<-rbind.fill(quantframe)

comm.write.table(meanC,"beta_quantiles.csv",row.names=FALSE)

finalize()

