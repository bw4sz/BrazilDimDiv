#Find unique richnesses in the matrix

#Require libraries
suppressMessages(require(data.table,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(picante,quietly=TRUE,warn.conflicts=FALSE))

##set path
droppath<-"/work/02443/bw4sz/GlobalMammals/"
#droppath<-"C:/Users/Caterina/Dropbox/R/1312_beta"
setwd(droppath)

########################################
###############Read in data on the node 0
########################################

##Read siteXspp table
siteXspp <- fread("Input/siteXspp1dgr.csv")

#make first identification rowname
siteXspp[,V1:=1:nrow(siteXspp)]

#match with phylogeny
tree<-read.tree("Input/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")
siteXspp<-siteXspp[,colnames(siteXspp) %in% c(tree$tip.label,"V1"),with=F]

#calculate richness
siteXspp$rich<-rowSums(siteXspp[,!colnames(siteXspp) %in% "V1",with=F])

#setkey to richness
setkey(siteXspp,rich)

#Remove lines with less than 2 species
siteXspp<-siteXspp[rich>1]

#Keep only original row (V1) and richness
siteXrich<-siteXspp[,colnames(siteXspp) %in% c("V1","rich"),with=F]

#write siteXrich table
write.csv(siteXrich,"Output/siteXrich.csv")

#remove duplicates
dt.unique<-unique(siteXrich)

#Get entire species list
splist<-colnames(siteXspp[, c("x","y","rich","V1"):=NULL,])

rm(siteXspp)

write.csv(dt.unique,"Input/UniquesiteXrich.csv",row.names=FALSE)