require(data.table)

#setwd
droppath<-"/home1/02443/bw4sz/GlobalMammals/"

#Betadiversity data
dat<-fread("Output/FinalData.txt")

print(dat)

#remove unwanted columns

dat[,To.x:=NULL]
dat[,From.x:=NULL]
dat[,To.y:=NULL]
dat[,From.y:=NULL]
dat[,combo:=NULL]

#xydist dist table

xydist<-data.table(read.table("Output/xydist.txt",row.names=1,header=TRUE))
print(xydist)


setnames(xydist,colnames(xydist),c("To.OriginalRow","From.OriginalRow","km"))

print(head(xydist))

setkey(dat,To.OriginalRow,From.OriginalRow)
setkey(xydist,To.OriginalRow,From.OriginalRow)

a<-merge(xydist,dat)

print(a)

setkey(xydist,From.OriginalRow,To.OriginalRow)

## reverse names
setnames(xydist,colnames(xydist),c("From.OriginalRow","To.OriginalRow","km"))

# set keys
setkey(dat,To.OriginalRow,From.OriginalRow)
setkey(xydist,To.OriginalRow,From.OriginalRow)

b<-merge(xydist,dat)

d<-rbind(a,b)

d<-unique(d)

print(d)

write.csv(d,"Output/BetaDist.csv",row.names=TRUE)


env<-data.table(read.table("Output/EnvData.txt",header=TRUE,row.names=NULL))

setnames(env,colnames(env),c("rowID","To.OriginalRow","From.OriginalRow","envdist"))

print(env)

setkey(dat,To.OriginalRow,From.OriginalRow)
setkey(env,To.OriginalRow,From.OriginalRow)

a<-merge(env,d)

print(a)

setkey(env,From.OriginalRow,To.OriginalRow)

## reverse names
setnames(env,colnames(env),c("rowID","From.OriginalRow","To.OriginalRow","envdist"))

# set keys

setkey(d,To.OriginalRow,From.OriginalRow)
setkey(env,To.OriginalRow,From.OriginalRow)

b<-merge(env,d)

BetaEnvDist<-rbind(a,b)

BetaEnvDist<-unique(BetaEnvDist)

print(BetaEnvDist)

write.table(BetaEnvDist,"Output/BetaDistEnv.txt",row.names=FALSE)
