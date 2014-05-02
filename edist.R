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



