require(data.table)

#setwd
droppath<-"/work/02443/bw4sz/GlobalMammals/"

# #Betadiversity data ##Cp from here to Line 146 it was deleted so I removed all the "#"
dat<-fread("Output/FinalData.txt")
setnames(dat,colnames(dat),c("To","From","Phylo","Tax","Trait"))

print("Beta")
print(dat)

###Step 1## Match unique ID to rowID and xytable
# ############################################################################
xytable<-fread("Output/xytable.csv")[,-1,with=F]

print(xytable)

# create all combinations of 2
# return as a data.table with these as columns `V1` and `V2`
# then count the numbers in each group

system.time(combs<-as.data.table(t(combn(xytable$V1,2, simplify = T))))

print(combs)

print(gc())

#Merge xy info, set keys
setkey(combs,V1)
setkey(xytable,V1)

#combine both ways
combsV1<-xytable[combs]
setkey(combsV1,V2)
combsV2<-xytable[combsV1]

#set names
setnames(combsV2,colnames(combsV2),c("From.OriginalRow","From.xcord","From.ycord","From","To.OriginalRow","To.xcord","To.ycord","To")) #CP inverted to and from because to is 1:n and from 1:n+1
#setnames(combsV2,colnames(combsV2),c("To.OriginalRow","To.xcord","To.ycord","To","From.OriginalRow","From.xcord","From.ycord","From"))

print("xytable")
print(combsV2)
rm(combs)

#########################################################################
#Merge with betadiversity calculations

#Bring in pairwise betadiveristy data
print(head(dat))

#set keys to match
setkey(dat,To,From)

#add all combinations that are 0 because To-From are the same
Todat<-dat[,"To",with=F]
Fromdat<-dat[,"From",with=F]
ToFrom<-unique(rbind(Todat,Fromdat,use.names=F))
suppdata<-cbind(ToFrom,ToFrom)
setnames(suppdata,1:2,c("To","From"))
betas0<-data.table(Phylo=rep(0,nrow(suppdata)),Tax=rep(0,nrow(suppdata)),Trait=rep(0,nrow(suppdata)))
suppdata2<-cbind(suppdata,betas0)
dat2<-rbind(dat,suppdata2)
rm(dat)
rm(betas0)
rm(suppdata)
rm(suppdata2)
rm(Fromdat)
rm(Todat)
rm(ToFrom)
gc()

#set to character, just help makes everyone match
combsV2[,To:=as.character(To)]
combsV2[,From:=as.character(From)]

dat2[,To:=as.character(To)]
dat2[,From:=as.character(From)]

#now set key
setkey(combsV2,To,From)

#first merge
mergeL<-merge(combsV2,dat2)

print("mergeL")

print(mergeL)

#Switch to and from column name 
setnames(dat2,colnames(dat2),c("From","To","Phylo","Tax","Trait"))
setkey(dat2,To,From)

#merge again
mergeR<-merge(combsV2,dat2)

print("mergeR")
print(mergeR)

print(gc())

rm(dat2)
rm(combsV2)
rm(combsV1)
rm(xytable)

tables()

print(gc())

datxy<-rbindlist(list(mergeR,mergeL)) #datxy has 72 662 658 rows but combsV2 has 75 614 253 rows
# it is just because some cells have been removed after creating xytable (cells having less than 2 species)

print("datxy")
print(datxy)

#remove unwanted columns
datxy[,To:=NULL]
datxy[,From:=NULL]

rm(mergeR)
rm(mergeL)
ls()

######################################################################
#Merge with xydistance
#xydist dist table

xydist<-fread("Output/xydist.txt")[,-1,with=F]

setnames(xydist,colnames(xydist),c("To.OriginalRow","From.OriginalRow","km"))

print(xydist)

setkey(datxy,To.OriginalRow,From.OriginalRow)
setkey(xydist,To.OriginalRow,From.OriginalRow)

a<-merge(xydist,datxy)

print(a)

#reverse names
setnames(xydist,colnames(xydist),c("From.OriginalRow","To.OriginalRow","km"))

#set keys
setkey(xydist,To.OriginalRow,From.OriginalRow)

b<-merge(xydist,datxy)

d<-rbindlist(list(a,b))

rm(a)
rm(b)
rm(datxy)
rm(xydist)

print(gc())

d<-unique(d)

print(d)

save.image("Output/BetaImage.RData")

############################################################################
#Merge with environmental betadiversity

load("Output/BetaImage.RData")


env<-fread("Output/EnvData.txt")[,-1,with=F]

setnames(env,colnames(env),c("To.OriginalRow","From.OriginalRow","envdist"))

print(env)

#set keys
setkey(d,To.OriginalRow,From.OriginalRow)
setkey(env,To.OriginalRow,From.OriginalRow)

#merge
a<-merge(env,d)

print(a)

setkey(env,From.OriginalRow,To.OriginalRow)

## reverse names
setnames(env,colnames(env),c("From.OriginalRow","To.OriginalRow","envdist"))

# set keys

setkey(d,To.OriginalRow,From.OriginalRow)
setkey(env,To.OriginalRow,From.OriginalRow)

b<-merge(env,d)

BetaEnvDist<-rbindlist(list(a,b))
rm(a)
rm(b)
rm(d)
rm(env)
print(gc())

BetaEnvDist<-unique(BetaEnvDist)

print(BetaEnvDist)

save.image("Output/BetaDistEnv.RData")

write.table(BetaEnvDist,"Output/BetaDistEnv.txt",row.names=FALSE)
