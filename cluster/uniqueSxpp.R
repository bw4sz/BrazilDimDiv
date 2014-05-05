#Find unique rows in the matrix

#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW


##################


#Require libraries, i hate the messages
suppressMessages(require(data.table,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(picante,quietly=TRUE,warn.conflicts=FALSE))

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

#Remove lines with less than 2 species
richness<-rowSums(dt.unique[,!colnames(dt.unique) %in% "id",with=F])
keep<-dt.unique[richness>1,"id",with=F]$id
dt.unique<-dt.unique[keep]

#Get entire species list
splist<-colnames(dt.unique)

#Read in phylogeny
tree<-read.tree("Input/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")

#remove species in the siteXspp that are not in phylogeny
dt.unique<-dt.unique[,colnames(dt.unique) %in% c(tree$tip.label,"id"),with=F]

print(dim(dt.unique))

print(rownames(dt.unique)[1:10])

print(paste("dimensions of the unique matrix:",dim(dt.unique)))

write.csv(dt.unique,"Input/UniquesiteXspp.csv",row.names=FALSE)