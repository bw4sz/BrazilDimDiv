require(data.table)
require(gRbase)

#setwd
droppath<-"/home1/02443/bw4sz/GlobalMammals/"

#Betadiversity data
dat<-fread("Output/FinalData.txt")

###Step 1## Match unique ID to rowID and xytable
############################################################################
xytable<-fread("Output/xytable.csv")[,-1,with=F]

# create all combinations of 2
# return as a data.table with these as columns `V1` and `V2`
# then count the numbers in each group

system.time(combs<-combn(xytable$V1,2,simplify = F))

combs<-do.call(rbind,combs)

combs<-data.table(combs)

#Merge xy info, set keys
setkey(combs,V1)
setkey(xytable,V1)

#combine both ways
combsV1<-xytable[combs]
setkey(combsV1,V2)
combsV2<-xytable[combsV1]

#set names
setnames(combsV2,colnames(combsV2),c("From.OriginalRow","From.xcord","From.ycord","From","To.OriginalRow","To.xcord","To.ycord","To")) #CP inverted to and from because to is 1:1233 and from 1:1234
#setnames(combsV2,colnames(combsV2),c("To.OriginalRow","To.xcord","To.ycord","To","From.OriginalRow","From.xcord","From.ycord","From"))

print(head(combsV2))

#########################################################################
#Merge with betadiversity calculations #CP these datasets should be merged by original row and not on to/from

#Bring in pairwise betadiveristy data
print(head(dat))

#set keys to match
#setkey(dat,To,From) #CP moved after set to character

#set to character, just help makes everyone match
combsV2[,To.OriginalRow:=as.character(To.OriginalRow)] #CP original rows
combsV2[,From.OriginalRow:=as.character(From.OriginalRow)] #CP original rows

dat[,To.OriginalRow:=as.character(To.OriginalRow)] #CP original rows
dat[,From.OriginalRow:=as.character(From.OriginalRow)] #CP original rows

#now set key
setkey(dat,To.OriginalRow,From.OriginalRow) #CP original rows
setkey(combsV2,To.OriginalRow,From.OriginalRow) #CP original rows

#first merge
mergeL<-merge(combsV2,dat)

#Switch to and from column name #CP is this needed? (and to goes from 1 to x and from goes from 2 to x+1)
setnames(combsV2,colnames(combsV2),c("To.OriginalRow","To.xcord","To.ycord","From","From.OriginalRow","From.xcord","From.ycord","To"))

#merge again
mergeR<-merge(combsV2,dat)

datxy<-rbind(mergeR,mergeL)

#remove unwanted columns
datxy[,To:=NULL]
datxy[,From:=NULL]

######################################################################
#Merge with xydistance
#xydist dist table

xydist<-fread("Output/xydist.txt")[,-1,with=F]
print(xydist)

setnames(xydist,colnames(xydist),c("To.OriginalRow","From.OriginalRow","km"))

print(head(xydist))

xydist[,To.OriginalRow:=as.character(To.OriginalRow)] # CP set to character 
xydist[,From.OriginalRow:=as.character(From.OriginalRow)] # CP set to character

setkey(datxy,To.OriginalRow,From.OriginalRow)
setkey(xydist,To.OriginalRow,From.OriginalRow)

a<-merge(xydist,datxy)

print(a)

## reverse names
setnames(xydist,colnames(xydist),c("From.OriginalRow","To.OriginalRow","km"))

# set keys
setkey(xydist,To.OriginalRow,From.OriginalRow)

b<-merge(xydist,datxy)

d<-rbind(a,b)

d<-unique(d)

print(d)

write.csv(d,"Output/BetaDist.csv",row.names=TRUE)
############################################################################
#Merge with environmental betadiversity

env<-data.table(read.table("Output/EnvData.txt",header=TRUE,row.names=NULL))

setnames(env,colnames(env),c("To.OriginalRow","From.OriginalRow","envdist")) # CP removed rowID
env[,To.OriginalRow:=as.character(To.OriginalRow)] # CP set to character 
env[,From.OriginalRow:=as.character(From.OriginalRow)] # CP set to character

print(env)

#set keys
setkey(d,To.OriginalRow,From.OriginalRow) #CP replaced dat with d
setkey(env,To.OriginalRow,From.OriginalRow)

#merge
a<-merge(env,d)

print(a)

#setkey(env,From.OriginalRow,To.OriginalRow) #CP delete

## reverse names
setnames(env,colnames(env),c("rowID","From.OriginalRow","To.OriginalRow","envdist"))

# set keys

setkey(d,To.OriginalRow,From.OriginalRow)
setkey(env,To.OriginalRow,From.OriginalRow)

b<-merge(env,d)

BetaEnvDist<-rbind(a,b)

BetaEnvDist<-unique(BetaEnvDist)

print(BetaEnvDist)

BetaEnvDist2<-BetaEnvDist[,c(1:10,13:15),,with=F] #CP keep only unique columns


write.table(BetaEnvDist,"Output/BetaDistEnv2.txt",row.names=FALSE) #CP name changed
