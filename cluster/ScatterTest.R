#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW


##################

require(pbdMPI,quietly=TRUE,warn.conflicts=FALSE)

init()

#Require libraries, i hate the messages
suppressMessages(require(reshape,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(picante,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(ggplot2,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(parallel,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(foreach,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(GGally,quietly=TRUE,warn.conflicts=FALSE))
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

if (comm.rank()==0){
siteXspp <- fread("Input/siteXspp1dgr.csv")    # as usual R read.table

#make v1 a key column
siteXspp[,V1:=1:nrow(siteXspp)]
setkey(siteXspp,"V1")

print(paste(comm.rank(),"loaded"))

#Read in phylogeny
tree<-read.tree("Input/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")

#bring in traits
traits <- read.table("Input/imputedmammals28apr14.txt",header=TRUE,row.names=1)

print(head(traits))
print(dim(traits))

#Read in cell by branch table made from source
load(file="Input/tcellbr.RData")
print("read in tcellbr.RData")

print("data loaded!")

#Remove xy data
siteXspp<-siteXspp[, c("x","y","rich"):=NULL,]

#Remove lines with less than 2 species
system.time(richness<-rowSums(siteXspp[,-1,with=F]))

keep<-which(richness > 1)

#Keep siteXspp columns
siteXspp<-siteXspp[keep,1:ncol(siteXspp),with=F]

#Get entire species list
splist<-colnames(siteXspp)

#remove species in the siteXspp that are not in phylogeny
siteXspp<-siteXspp[,c(colnames(siteXspp)[colnames(siteXspp) %in% tree$tip.label],key(siteXspp)),with=F]
print(dim(siteXspp))
print(dim(siteXspp))
print(paste("There are in siteXspp NAs:",sum(is.na(siteXspp))))

#subtest
system.time(comm<-siteXspp[V1 <= 500])

print(paste("There are in comm NAs:",sum(is.na(comm))))

#Create all pairwise combinations of siteXspp
z<-combn(comm$V1,2)

#index by key
print(paste("Total number of iterations:",ncol(z)))

#split rows into indices
IndexFunction<-splitIndices(ncol(z),comm.size())

#print(paste("Length of IndexFunction is:",length(IndexFunction)))

  ###Send the subset matrix
  toScatterMatrix<-lapply(IndexFunction,function(y){
    Index_Space<-z[,y]
    rowsTocall<-unique(as.vector(Index_Space))
    comm.df<-data.frame(comm[V1 %in% rowsTocall])
    comm.df<-comm.df[,which(!apply(comm.df,2,sum)==0)]
    rownames(comm.df)<-comm.df$V1
    comm.df<-comm.df[,!colnames(comm.df) %in% "V1"]
  })
  
  ###Send the subset 
  toScatterIndex<-lapply(IndexFunction,function(y){
    Index_Space<-z[,y]
  })
  
  ##subset of the tcellbr
  toScatterTcell<-lapply(IndexFunction,function(y){
    Index_Space<-z[,y]
    rowsTocall<-unique(as.vector(Index_Space))
    a<-tcellbr[rownames(tcellbr) %in% rowsTocall,]
  })

##subset of the tcellbr rownames 
toScatterTcellrownames<-lapply(IndexFunction,function(y){
  Index_Space<-z[,y]
  rowsTocall<-unique(as.vector(Index_Space))
  rownames(tcellbr)[rownames(tcellbr) %in% rowsTocall]
})


} else  {
toScatterTcellrownames<-NULL
traits<-NULL
toScatterMatrix<-NULL
toScatterIndex<-NULL
toScatterTcell<-NULL
}


#get data for that core
dat<-scatter(toScatterMatrix,rank.source=0)

#get index for that core
Indat<-scatter(toScatterIndex,rank.source=0)

#get tcell for that core
Tcell<-scatter(toScatterTcell,rank.source=0)

Tcellnames<-scatter(toScatterTcellrownames,rank.source=0)

#Cast trait matrix to everyone
traits<-bcast(traits)

comm.print(Indat[,1:5],all.rank=TRUE)

comm.print(Tcell[1:5,1:5],all.rank=TRUE)


comm.print(paste("Total nodes is:", comm.size()))

timeF<-system.time(beta_out<-betaPar.scatter(toScatterMatrix=dat,toScatterIndex=Indat,tcellbr=Tcell,rown=Tcellnames))

#Return timing argument to console
print(timeF)

#try writing from all
comm.write.table(beta_out,"beta_out.csv")

finalize()

#Data Generation Complete
##########################################################################################
##########################################################################################


finalize()