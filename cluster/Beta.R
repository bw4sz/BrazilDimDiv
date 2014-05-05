#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW


##################

require(pbdMPI,quietly=TRUE,warn.conflicts=FALSE)

init()

#create memory profile

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
siteXspp <- fread("Input/UniquesiteXspp.csv")    # as usual R read.table

print("siteXspp table:")
print(siteXspp[1:5,1:5,with=F])

#make v1 a key column
siteXspp[,V1:=1:nrow(siteXspp)]
setkey(siteXspp,V1)

print("Tables in memory:")
tables()

#Read in phylogeny
tree<-read.tree("Input/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")

#bring in traits
traits <- read.table("Input/imputedmammals28apr14.txt",header=TRUE,row.names=1)

#Read in cell by branch table made from source
load(file="Input/tcellbr.RData")

#subtest
comm<-siteXspp[V1 < 10]

rm(siteXspp)

#Create all pairwise combinations of siteXspp
z<-combn(comm$V1,2)

#index by key
print(paste("Total number of iterations:",ncol(z)))

#split rows into indices
IndexFunction<-splitIndices(ncol(z),comm.size())

print(paste("Length of IndexFunction is:",length(IndexFunction)))

  ###Send the subset matrix
  toScatterMatrix<-lapply(IndexFunction,function(y){
    Index_Space<-z[,y]
    rowsTocall<-unique(as.vector(Index_Space))
    comm.df<-data.frame(comm[V1 %in% rowsTocall])
    comm.df<-comm.df[,which(!apply(comm.df,2,sum)==0)]
    rownames(comm.df)<-comm.df$V1
    comm.df<-comm.df[,!colnames(comm.df) %in% "V1"]
  })
  
rm(comm)

print("toScatterMatrix")
  ###Send the subset 
  toScatterIndex<-lapply(IndexFunction,function(y){
    Index_Space<-z[,y]
  })
print("toScatterIndex")

  
  ##subset of the tcellbr
  toScatterTcell<-lapply(IndexFunction,function(y){
    Index_Space<-z[,y]
    rowsTocall<-unique(as.vector(Index_Space))
    a<-tcellbr[rownames(tcellbr) %in% rowsTocall,]
    b<-a[,which(!apply(a,2,sum)==0)]
return(b)
  })


print("toScatterTcell")

##subset of the tcellbr rownames 
toScatterTcellrownames<-lapply(IndexFunction,function(y){
  Index_Space<-z[,y]
  rowsTocall<-unique(as.vector(Index_Space))
  a<-tcellbr[rownames(tcellbr) %in% rowsTocall,]
    b<-a[,which(!apply(a,2,sum)==0)]
rownames(b)
})


print("toScatterTcellnames")

rm(tcellbr)


toScatterTrait<-lapply(toScatterMatrix,function(y){
traits[rownames(traits) %in% colnames(y),] 
})

print("toScatterTraits")

rm(traits)


} else  {
toScatterTcellrownames<-NULL
toScatterMatrix<-NULL
toScatterIndex<-NULL
toScatterTcell<-NULL
toScatterTrait<-NULL
}


#get data for that core
time_scatter<-system.time(dat<-scatter(toScatterMatrix,rank.source=0))

print(paste("Time to scatter matrix:",time_scatter[3]))
print(gc())

comm.print("scattered matrix")

#get index for that core
Indat<-scatter(toScatterIndex,rank.source=0)

comm.print("scattered index")

#get tcell for that core
Tcell<-scatter(toScatterTcell,rank.source=0)

comm.print("scattered tcell")

Tcellnames<-scatter(toScatterTcellrownames,rank.source=0)
comm.print("scattered tcellnames")

#Cast trait matrix to everyone
traitn<-scatter(toScatterTrait,rank.source=0)

comm.print("scattered trait")

comm.print(paste("Total nodes is:", comm.size()))

timeF<-system.time(beta_out<-betaPar.scatter(toScatterMatrix=dat,toScatterIndex=Indat,tcellbr=Tcell,rown=Tcellnames,traitn))

#Return timing argument to console
print(timeF)

#try writing from all
comm.write.table(beta_out,"beta_out.csv",row.names=FALSE)

finalize()

#Data Generation Complete
##########################################################################################
##########################################################################################


finalize()