library(pbdMPI,quietly=TRUE,warn.conflicts=FALSE,verbose=FALSE)

init()

#set position within the loop
comm.print(paste("Total nodes is:", comm.size()))

#set rank
.rank<-comm.rank()+1

#load libraries
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
splist<-colnames(siteXspp[, c("rich","V1"):=NULL,])

rm(siteXspp)

##We need a null model from every richness type to every other richness type. Including to itself.
ident<-expand.grid(dt.unique$rich,dt.unique$rich)

#Chunk the identity table based on the number of ranks
IndexFunction<-splitIndices(nrow(ident),comm.size())

#Which rank am I?
myIndexes<-IndexFunction[[.rank]]

##############Make replicates of your rank###############
#For each row of null assemblages, make 100 samples

#make a list
nullIndex<-list()

#number of replicates per assemblage type
iter<-5

#make a holder
quantile_output<-list()

#For each set of null assembles sample assemblage and run betascatter
for (row in 1:length(myIndexes)){
  
  outA<-list()
  for (g in 1:iter){
    outA[[g]]<-sampleS(myIndexes[row],ident)  
  }
  
  comm.df<-rbind.fill(outA)

  #Fill the absences to 0
  comm.df[is.na(comm.df)]<-0

  #Calculate phylogenetic species lists
  sp.list.phylo<-apply(comm.df,1,function(x){ 
    names(x[which(x==1)])
  })

#remove any species for which we have no trait information
  sp.list.trait<-apply(comm.df,1,function(x){ 
    trait_sp<-names(x[which(x==1)])
    trait_sp[trait_sp %in% colnames(traitdistance)]
  })


#Which rows to compare, we have randomly sampled, so just go in order, 
Index_Space<-t(data.frame(To=seq(1,nrow(comm.df),2),From=seq(2,nrow(comm.df),2)))

#needed be character names, not numbers
names(sp.list.trait)<-as.character(1:length(sp.list.trait))
names(sp.list.phylo)<-as.character(1:length(sp.list.phylo))

beta_out<-betaPar.scatter(toScatterIndex = Index_Space,coph=coph,traitdist=traitdistance,sp.list.phylo = sp.list.phylo,sp.list.trait = sp.list.trait)

#get quantiles for each type
quants<-as.data.frame(apply(beta_out[,-c(1,2)],2,function(x){
  c(Lower=quantile(x,0.025),Upper=quantile(x,0.975))
}))

quants$Quantile<-c("Lower","Upper")
rownames(quants)<-NULL

#What richness type was this
quants$RichnessA<-ident[myIndexes[row],1]
quants$RichnessB<-ident[myIndexes[row],2]

quantile_output[[row]]<-quants
}

#Bind all the outputs together
quantile_all<-rbind.fill(quantile_output)

#checkpoint, write data if fails
save.image(paste("Beta/",.rank,"Randomization.RData",sep=""))

#Return timing argument to console
#print(timeF)

#try writing from all
comm.write.table(beta_out,"Output/Randomization.txt",row.names=F,append=T,col.names = FALSE)


#To delete, but if you want to visualize it do this
dat<-melt(quantile_all,id.vars=c("RichnessA","RichnessB","Quantile"))
dat<-dcast(dat,...~Quantile)

ggplot(dat,aes(x=paste(RichnessA,"_",RichnessB),col=variable)) + geom_linerange(aes(ymin=Lower,ymax=Upper)) + labs(x="Richness") + facet_wrap(~variable,scale="free")

#remove data file
#file.remove(fil)

finalize()

###############################
#Null Generation Complete
###############################
