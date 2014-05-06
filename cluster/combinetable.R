#require packages
require(data.table)

#set working directory
droppath<-"/home1/02443/bw4sz/GlobalMammals/"
setwd(droppath)

#read in xy data with original rows
xytab<-fread("xytable.csv",nrows=1000)[,-1,with=F]
setnames(xytab,colnames(xytab),c("V1","x","y","id2"))

# create all combinations of 2
# return as a data.table with these as columns `V1` and `V2`
# then count the numbers in each group
combs<-xytab[, rbindlist(combn(V1,2, 
                      FUN = function(x) as.data.table(as.list(x)), simplify = F))]
setkey(combs,V1)
setkey(xytab,V1)

#combine both ways
combsV1<-xytab[combs]
setkey(combsV1,V2)
combsV2<-xytab[combsV1]

#set names
setnames(combsV2,colnames(combsV2),c("To.OriginalRow","To.x","To.y","To","From.OriginalRow","From.x","From.y","From"))

#Bring in pairwise betadiveristy data
beta_out<-fread("beta_out.csv")
print(head(beta_out))

#set keys to match
setkey(beta_out,To,From)

#set to character
combsV2[,To:=as.character(To)]
combsV2[,From:=as.character(From)]

setkey(combsV2,To,From)

head(mergeT<-merge(combsV2,beta_out))

print(dim(mergeT))

write.csv(Frommerge,"FinalData.csv")
