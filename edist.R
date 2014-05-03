##Our goal is to compute efficient pairwise binary calculations.

#Set seed
set.seed(123)

#Create a data table with a key column and a bunch a 01 columns (sp*)
dt<-data.table(V1=1:10,spA=rbinom(10,1,.5),spB=rbinom(10,1,.5))

##We can notice a few things a
setkey(dt,V1)

#Which are the unique combinations of columns, besides the key
#This is a bit ugly, suggestions welcome
dt.unique<-subset(dt,!duplicated(dt[,!colnames(dt) %in% "V1",with=F]))

#compute distance matrix on unique subset
dist.unique<-as.matrix(dist(dt.unique[,!colnames(dt.unique) %in% "V1",with=F]))

#melt into a data.frame
distmelt<-melt(dist.unique)

#name, key and make a data.table
colnames(distmelt)<-c("KeyTo","KeyFrom","dist")

dtm<-data.table(distmelt)
#Create combination key
dtm[,Combo:=paste(KeyTo,KeyFrom)]
setkey(dtm,KeyTo,KeyFrom)

#create combination table to hold output
combs<-t(combn(dt$V1,2))
colnames(combs)<-c("KeyTo","KeyFrom")

#make into data.table
dt.out<-data.table(combs)
setkey(dt.out,KeyTo,KeyFrom)

#merge combinations with distance table
dt.merge<-dt.out[dtm]

#remove the diagonal these were not calcuated, but created with melting the table, ideas appreciated
dt.merge[!KeyTo==KeyFrom]

#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW


##################


#Require libraries, i hate the messages
suppressMessages(require(data.table,quietly=TRUE,warn.conflicts=FALSE))

#Everyone say hello
#comm.print(comm.rank(), all.rank = TRUE)

##If running locally set droppath

droppath<-"C:/Users/Ben/Documents/BrazilDimDiv/cluster"

setwd(droppath)

###Define Source Functions, does this need to be run and distributed to all nodes, can they source simultaneously

suppressMessages(source("Input/BrazilSourceFunctions.R"))


########################################
###############Read in data on the node 0
########################################

siteXspp <- fread("Input/siteXspp1dgr.csv")    # as usual R read.table

#make v1 a key column
siteXspp[,V1:=1:nrow(siteXspp)]
setkey(siteXspp,"V1")

#Remove xy data
siteXspp<-siteXspp[, c("x","y","rich"):=NULL,]

#Remove lines with less than 2 species
system.time(richness<-rowSums(siteXspp[,-1,with=F]))

keep<-which(richness > 1)

#Keep siteXspp columns
siteXspp<-siteXspp[keep,1:ncol(siteXspp),with=F]

#subtest
system.time(comm<-siteXspp[V1 < 3000])

print(paste("There are in comm NAs:",sum(is.na(comm))))

#Create all pairwise combinations of siteXspp
z<-combn(comm$V1,2)

#index by key
print(paste("Total number of iterations:",ncol(z)))

#split rows into indices
IndexFunction<-splitIndices(ncol(z),comm.size())

print(paste("Length of IndexFunction is:",length(IndexFunction)))

###Send the subset matrix
toScatterMatrix<-lapply(IndexFunction,function(y){
  Index_Space<-z[,y]
  rowsTocall<-unique(as.vector(Index_Space))
  comm.df<-data.frame(comm[V1 %in% rowsTocall])
  comm.df<-comm.df[,which(!apply(comm.df,2,sum)==0)]
  rownames(comm.df)<-comm.df$V1
  comm.df<-comm.df[,!colnames(comm.df) %in% "V1"]
})

rm(comm)

gc()

print("toScatterMatrix")
###Send the subset 
toScatterIndex<-lapply(IndexFunction,function(y){
  Index_Space<-z[,y]
})
print("toScatterIndex")


##subset of the tcellbr
toScatterTcell<-lapply(IndexFunction,function(y){
  Index_Space<-z[,y]
  rowsTocall<-unique(as.vector(Index_Space))
  a<-tcellbr[rownames(tcellbr) %in% rowsTocall,]
  b<-a[,which(!apply(a,2,sum)==0)]
})


print("toScatterTcell")

##subset of the tcellbr rownames 
toScatterTcellrownames<-lapply(IndexFunction,function(y){
  Index_Space<-z[,y]
  rowsTocall<-unique(as.vector(Index_Space))
  a<-tcellbr[rownames(tcellbr) %in% rowsTocall,]
  b<-a[,which(!apply(a,2,sum)==0)]
  rownames(b)
})


print("toScatterTcellnames")

rm(tcellbr)

gc()


toScatterTrait<-lapply(toScatterMatrix,function(y){
  traits[rownames(traits) %in% colnames(y),] 
})

print("toScatterTraits")

rm(traits)


} else  {
  toScatterTcellrownames<-NULL
  toScatterMatrix<-NULL
  toScatterIndex<-NULL
  toScatterTcell<-NULL
  toScatterTrait<-NULL
}



#get data for that core
dat<-scatter(toScatterMatrix,rank.source=0)

comm.print("scatter matrix")

#get index for that core
Indat<-scatter(toScatterIndex,rank.source=0)

comm.print("scatter index")

#get tcell for that core
Tcell<-scatter(toScatterTcell,rank.source=0)

comm.print("scatter tcell")

Tcellnames<-scatter(toScatterTcellrownames,rank.source=0)
comm.print("scatter tcellnames")

#Cast trait matrix to everyone
traitn<-scatter(toScatterTrait,rank.source=0)

comm.print("scatter trait")


comm.print(paste("Total nodes is:", comm.size()))

timeF<-system.time(beta_out<-betaPar.scatter(toScatterMatrix=dat,toScatterIndex=Indat,tcellbr=Tcell,rown=Tcellnames,traitn))

#Return timing argument to console
print(timeF)

#try writing from all
comm.write.table(beta_out,"beta_out.csv")

finalize()

#Data Generation Complete
##########################################################################################
##########################################################################################


finalize()


