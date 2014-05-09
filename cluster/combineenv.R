require(data.table)

#setwd
setwd("cluster")

#Betadiversity data
dat<-fread("Output/FinalData.csv",verbose=TRUE)

print(head(dat))

setkey(dat,To.OriginalRow,From.OriginalRow)
#env dist table

env<-fread("Output/EnvData.csv",verbose=TRUE)[,-1,with=F]

setnames(env,colnames(env),c("To.OriginalRow","From.OriginalRow","Env.dist"))

setkey(env,To.OriginalRow,From.OriginalRow)

#merge betadiversity and env diversity data.

#make sure both tables have the same keys
tables()

head(mergeT<-merge(dat,env))

