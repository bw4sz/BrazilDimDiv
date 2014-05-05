###Read 

require(data.table)

droppath<-"/home1/02443/bw4sz/GlobalMammals/"

setwd(droppath)

xytab<-fread("xytable.csv",nrows=10)[,-1,with=F]

setnames(xytab,colnames(xytab),c("V1","x","y","id2"))

xytab[,id:=as.character(xytab$id)]
setkey(xytab,id)

# create all combinations of 2
# return as a data.table with these as columns `V1` and `V2`
# then count the numbers in each group
xytab[, rbindlist(combn(V1,2, 
                      FUN = function(x) as.data.table(as.list(x)), simplify = F))]


comb.table<-data.table(xytab[V1 %in% combs[,1]],xytab[V1 %in% combs[,2]])
setnames(comb.table,colnames(comb.table),c("To.Original",'To.x',"To.y","To.id","From.Original",'From.x',"From.y","From.id"))

beta_out<-fread("beta_out.csv")

print(head(beta_out))

#set keys
setkey(beta_out,To)

#merge with original xytab
head(Tomerge<-beta_out[xytab])

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
