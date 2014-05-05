###Read 

require(data.table)

setwd("C:/Users/Ben/Documents/BrazilDimDiv/cluster")

xytab<-fread("xytable.csv",nrows=100)[,-1,with=F]

setnames(xytab,colnames(xytab),c("V1","x","y","id"))

xytab[,id:=as.character(xytab$id)]
setkey(xytab,id)

beta_out<-fread("C:/Users/Ben/Dropbox/Dimensions Data/beta_out.csv",nrows=200)[,-1,with=F]

print(head(beta_out))
setnames(beta_out,colnames(beta_out),c("To","From","Tax","Phylo",'Trait'))

#set keys
setkey(beta_out,To)

#merge with original xytab
head(Tomerge<-xytab[beta_out,])

#rename the columns
setnames(Tomerge,colnames(Tomerge),c("To.id","To.OriginalRow",'To.x',"To.y","From.id","Tax","Phylo","Trait"))

#repeat with from
setkey(Tomerge,From.id)

head(Frommerge<-Tomerge[xytab,])
setnames(Frommmerge,colnames(Frommerge),c("To.id","To.OriginalRow",'To.x',"To.y","From.id","Tax","Phylo","Trait","From.OriginalRow","From.x","From.y"))

Frommerge[complete.cases(Frommerge),]

write.csv("FinalData.csv")
