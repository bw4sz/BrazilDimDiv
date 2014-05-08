require(data.table)

#setwd
setwd("cluster")

#Betadiversity data
dat<-fread("FinalData.csv")

print(head(dat))

#env dist table

env<-read.table("EnvData.csv",nrows=1000)

env.d<-data.table(env)

setnames(env.d,colnames(env.d),c("To","From","Env.dist"))

setkey(env.d,To,From)

fread("Input/")