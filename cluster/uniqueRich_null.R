library(pbdMPI,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)

init()

#set position within the loop
comm.print(paste("Total nodes is:", comm.size()))

#set rank
.rank<-comm.rank()+1

#library libraries
suppressMessages(library(reshape2,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(picante,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(ggplot2,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(parallel,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(foreach,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(GGally,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(data.table,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(betapart,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(library(plyr,quietly=TRUE,warn.conflicts=FALSE))

#Everyone say hello
comm.print(comm.rank(), all.rank = TRUE)

##set path
#droppath<-"/work/02443/bw4sz/GlobalMammals/"
droppath<-"C:/Users/Ben/Documents/BrazilDimDiv/cluster/"
setwd(droppath)

###Define Source Functions
suppressMessages(source("Input/BrazilSourceFunctions.R"))

###Read in data

#Bring in relatedness matrix
coph<-fread("Input/cophenetic.csv")
setkey(coph,V1)
setcolorder(coph,c("V1",coph$V1))

#bring in trait distance matrix
traitdistance<-fread("Input/traitdistanceLog.csv")
setkey(traitdistance,V1)

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

#remove duplicates
dt.unique<-unique(siteXrich)

#Get entire species list
splist<-colnames(siteXspp[, c("x","y","rich","V1"):=NULL,])

rm(siteXspp)

##We need a null model from every richness type to every other richness type. Including to itself.
ident<-expand.grid(dt.unique$rich,dt.unique$rich)

#Chunk the identity table based on the number of ranks
IndexFunction<-splitIndices(nrow(ident),50)

#Which rank am I?
myIndexes<-IndexFunction[[.rank]]

##############Make replicates of your rank###############
#For each row of null assemblages, make 100 samples
outA<-list()
for (g in 1:100){
  outA[[g]]<-sampleS(myIndexes[1],ident)  
}

comm_df<-rbind.fill(outA)
#Fill the absences to 0
comm_df[is.na(comm_df)]<-0

#This would need to be done for every position in myIndexes, and they use betascatter to walk through the combinations,
#I have to stop here for the moment.

#Calculate phylogenetic species lists
sp.list.phylo<-apply(comm.df[,!colnames(comm.df) %in% "id",with=F],1,function(x){ 
  names(x[which(x==1)])
})

names(sp.list.phylo)<-comm.df$id

col_keep<-colnames(comm.df)[colnames(comm.df) %in% colnames(traitm)]
comm.df.trait<-comm.df[,c(col_keep,"id"),with=F]

sp.list.trait<-apply(comm.df.trait[,!colnames(comm.df.trait) %in% "id",with=F],1,function(x){ 
  names(x[which(x==1)])
})

names(sp.list.trait)<-comm.df$id

system.time(beta_out<-betaPar.scatter(toScatterIndex = Index_Space,coph=cophm,traitdist=traitm,sp.list.phylo = sp.list.phylo,sp.list.trait = sp.list.trait))

print(head(beta_out))

#checkpoint, write data if fails
save.image(paste("Beta/",.rank,"Randomization.RData",sep=""))

#Return timing argument to console
#print(timeF)

#try writing from all
comm.write.table(beta_out,"Output/Randomization.txt",row.names=F,append=T,col.names = FALSE)

#remove data file
#file.remove(fil)

finalize()

###############################
#Data Generation Complete
###############################

#now we run that pair 
out<-beta_all(coph=coph,traitdist=traitdistance,sp.list.phylo = splist,sp.list.trait = sp.list.trait,ids=x)
