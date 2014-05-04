#Find unique rows in the matrix

#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW


##################


#Require libraries, i hate the messages
suppressMessages(require(data.table,quietly=TRUE,warn.conflicts=FALSE))

#Everyone say hello
#comm.print(comm.rank(), all.rank = TRUE)

##If running locally set droppath

droppath<-"/home1/02443/bw4sz/GlobalMammals/"

setwd(droppath)

###Define Source Functions, does this need to be run and distributed to all nodes, can they source simultaneously

suppressMessages(source("Input/BrazilSourceFunctions.R"))


########################################
###############Read in data on the node 0
########################################
siteXspp <- fread("Input/siteXspp1dgr.csv")    # as usual R read.table

print(dim(siteXspp))

#make first identification rowname
siteXspp[,V1:=1:nrow(siteXspp)]

#setkey to all species
setkeyv(siteXspp,colnames(siteXspp)[!colnames(siteXspp) %in%  c("x","y","V1","rich")])

#make group counter
siteXspp[,id:=.GRP,by=key(siteXspp)]

#make xy v0 table
xytable<-siteXspp[,colnames(siteXspp) %in% c("x","y","V1","id"),with=F]

#write xytable
write.csv(xytable,"xytable.csv")

#Remove xy data
siteXspp<-siteXspp[, c("x","y","rich","V1"):=NULL,]

dt.unique<-subset(siteXspp,!duplicated(siteXspp))

write.csv(dt.unique,"Input/UniquesiteXspp.csv")