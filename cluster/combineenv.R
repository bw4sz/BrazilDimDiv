require(data.table)

#setwd
setwd("cluster")

#Betadiversity data
dat<-fread("Output/FinalData.csv",verbose=TRUE)

print(head(dat))

#env dist table

env<-fread("Output/EnvData.csv",verbose=TRUE)

rm(env)

gc()

setnames(env.d,colnames(env.d),c("To","From","Env.dist"))

setkey(env.d,To,From)

