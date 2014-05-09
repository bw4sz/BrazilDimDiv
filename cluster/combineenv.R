require(data.table)

#setwd
droppath<-"/home1/02443/bw4sz/GlobalMammals/"

#Betadiversity data
dat<-fread("Output/FinalData.csv",verbose=TRUE)

print(head(dat))

#remove unwanted columns

dat[,To.x:=NULL]
dat[,From.x:=NULL]
dat[,To.y:=NULL]
dat[,From.y:=NULL]
dat[,combo:=NULL]

#xydist dist table

xydist<-data.table(read.csv("Output/xydist.csv",row.names=1))

setnames(xydist,colnames(xydist),c("To.OriginalRow","From.OriginalRow","km"))

print(head(xydist))

setkey(dat,To.OriginalRow,From.OriginalRow)
setkey(xydist,To.OriginalRow,From.OriginalRow)

a<-merge(xydist,dat)

print(a)

setkey(xydist,From.OriginalRow,To.OriginalRow)

b<-merge(xydist,dat)

d<-rbind(a,b)

d<-unique(d)

print(d)

write.csv(d,"BetaDist.csv",row.names=TRUE)
