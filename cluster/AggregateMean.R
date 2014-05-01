#Compute average betadiversity for each cell, both in terms of mean value and quantiles for each pair. 

require(pbdMPI,quietly=TRUE,warn.conflicts=FALSE)

init()
require(parallel,quietly=TRUE,warn.conflicts=FALSE)
require(reshape,quietly=TRUE,warn.conflicts=FALSE)

#read in data

droppath<-"/home1/02443/bw4sz/GlobalMammals/"

#set rank number and size.
.rank<-comm.rank()
.size<-comm.size()

#read in data on rank 0

if(.rank==0){

#read in beta file from the cluster
dat5<-read.table(paste(droppath,"beta_out.csv",sep="/"),nrows=5,row.names=NULL)
classes <- sapply(dat5, class)

print(classes)

dat<-read.table(paste(droppath,"beta_out.csv",sep="/"),row.names=NULL)

print(head(dat))

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
comm.index<-indexS[[.rank+1]]

#loop through each in the index


meanC<-lapply(comm.index,function(x){
  
  #find cell value
  target<-cells[x]
  
  #find mean
  cellX<-dat[dat$To %in% target | dat$From %in% target,]
  meanCell<-t(sapply(cellX[,c("Sorenson","BetaSim","MNTD")],mean,na.rm=TRUE))
  return(data.frame(cbind(Cell=target,meanCell)))
})

meanCelltable<-rbind.fill(meanC)

comm.write.table(meanCelltable,"meanCelltable.txt")

#########################################
#Compute standard deviation of each cell
#########################################

meanC<-lapply(comm.index,function(x){
  
  #find cell value
  target<-cells[x]
  
  #find mean
  cellX<-dat[dat$To %in% target | dat$From %in% target,]
  meanCell<-t(sapply(cellX[,c("Sorenson","BetaSim","MNTD")],sd,na.rm=TRUE))
  return(data.frame(cbind(Cell=target,meanCell)))
})

meanCelltable<-rbind.fill(meanC)

comm.write.table(meanCelltable,"meanCelltable.txt")


##############################################
#Compute quantile for each combination of cells
##############################################

#for each row in data create a list position

indexAll<-splitIndices(nrow(dat),.size)

#have the cluster pull in based on rank number
comm.indexAll<-indexAll[[.rank+1]]

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

comm.write.table(beta_quantiles,"beta_quantiles.txt",row.names=FALSE)

finalize()

