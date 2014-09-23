suppressMessages(require(data.table))
suppressMessages(require(picante))
require(reshape2)


###############
#setwd on cluster
droppath<-"/work/02443/bw4sz/GlobalMammals/"

##If running locally set dropbox path
#droppath<-"C:/Users/Jorge/Documents/BrazilDimDiv/cluster"
setwd(droppath)

  env<-fread("Input/siteXsppXenv_1dgr.csv")
  xytable<-fread("xytable.csv")[,-1,with=F]

env.table<-env[,colnames(env) %in% c("x","y",paste("bio",1:19,sep="")),with=F]
setkey(env.table,x,y)
setkey(xytable,x,y)
env.table.id<-merge(xytable,env.table)
setkey(env.table.id,V1) 

env.pca<-prcomp(env.table.id[,-c(1:4),with=F])$x


xymat <- as.matrix(dist(env.pca))

print(dim(xymat))

colnames(xymat)<-env.table.id$V1
rownames(xymat)<-env.table.id$V1

#turn off diag and upper tri
diag(xymat)<-NA

xymat[upper.tri(xymat)]<-NA

xymelt<-melt(xymat)

xyout<-xymelt[!is.na(xymelt$value),]

#name the columns
colnames(xyout)<-c("To.OriginalRow","From.OriginalRow","envdist")

print(dim(xyout))

print(head(xyout))

write.table(xyout,"Output/EnvData.txt")