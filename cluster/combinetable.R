###Read 

require(data.table)

droppath<-"/home1/02443/bw4sz/GlobalMammals/"

setwd(droppath)

xytab<-fread("xytable.csv")[,-1,with=F]

setnames(xytab,colnames(xytab),c("V1","x","y","id"))

xytab[,id:=as.character(xytab$id)]
setkey(xytab,id)

beta_out<-fread("beta_out.csv")

print(head(beta_out))

#set keys
setkey(beta_out,To)

#merge with original xytab
head(Tomerge<-xytab[beta_out,])

#rename the columns
setnames(Tomerge,colnames(Tomerge),c("To.id","To.OriginalRow",'To.x',"To.y","From.id","BetaSim","Sorenson","MNTD"))

#repeat with from
setkey(Tomerge,From.id)

#buggy, repeat set to character
xytab[,id:=as.character(xytab$id)]
setkey(xytab,id)

head(Frommerge<-xytab[Tomerge,])
setnames(Frommerge,colnames(Frommerge),c("From.id","From.OriginalRow",'From.x',"From.y","To.id","To.OriginalRow","To.x","To.y","BetaSim","Sorenson","MNTD"))

Frommerge<-Frommerge[complete.cases(Frommerge),]


print(dim(Frommerge))

write.csv(Frommerge,"FinalData.csv")