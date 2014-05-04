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
#make v1 a key column
siteXspp[,V0:=1:nrow(siteXspp)]
setkey(siteXspp,"V0")

#Remove xy data
siteXspp<-siteXspp[, c("x","y","rich","V1"):=NULL,]

#Remove lines with less than 2 species
system.time(richness<-rowSums(siteXspp[,-1,with=F]))

keep<-which(richness > 1)

#Keep siteXspp columns
siteXspp<-siteXspp[keep,1:ncol(siteXspp),with=F]

#subtest
system.time(comm<-siteXspp)

print(paste("There are in comm NAs:",sum(is.na(comm))))
dt.unique<-subset(comm,!duplicated(comm[,!colnames(comm) %in% "V0",with=F]))

#make new column for new key
dt.unique[,V1:=1:nrow(dt.unique)]

print(dim(dt.unique))

write.csv(dt.unique,"Input/UniquesiteXspp.csv")
