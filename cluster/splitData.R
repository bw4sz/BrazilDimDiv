#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals
#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW


##################

library(pbdMPI,quietly=TRUE,warn.conflicts=FALSE)

init()

if(comm.rank()==0){
#Require libraries, i hate the messages
suppressMessages(require(reshape,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(picante,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(ggplot2,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(parallel,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(foreach,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(GGally,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(data.table,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(betapart,quietly=TRUE,warn.conflicts=FALSE))

#Everyone say hello
#comm.print(comm.rank(), all.rank = TRUE)

##If running locally set droppath

droppath<-"/work/02443/bw4sz/GlobalMammals/"

setwd(droppath)

###Define Source Functions, does this need to be run and distributed to all nodes, can they source simultaneously
comm <- fread("~/Input/UniquesiteXspp.csv")    # as usual R read.table

head(comm)

#read in xy data with original rows
xytab<-fread("Output/xytable.csv")[,-1,with=F]
setnames(xytab,colnames(xytab),c("V1","x","y","id"))

print("siteXspp table:")
print(comm[1:5,1:5,with=F])

#make v1 a key column
setkey(comm,id)

print("Tables in memory:")
tables()

#Bring in relatedness matrix
coph<-read.csv("Input/cophenetic.csv",row.names=1)

#bring in traits
traits <- read.table("Input/imputedmammals28apr14.txt",header=TRUE,row.names=1)

#bring in trait distance matrix
traitdistance<-read.csv("Input/traitdistance.csv",row.names=1)

#Create all pairwise combinations of siteXspp
z<-combn(comm$id,2)

#index by key
print(paste("Total number of iterations:",ncol(z)))

#split rows into indices
IndexFunction<-splitIndices(ncol(z),comm.size())
print(paste("Length of IndexFunction is:",length(IndexFunction)))

###Send the subset matrix

for (x in 1:length(IndexFunction)){
  print(x)
  #Index of rows to call.
  Index_Space<-z[,IndexFunction[[x]]]
  
  #row positions
  rowsTocall<-unique(as.vector(Index_Space))
  comm.df<-data.frame(comm[id %in% rowsTocall])
  comm.df<-comm.df[,which(!apply(comm.df,2,sum)==0)]
  rownames(comm.df)<-comm.df$id
  comm.df<-comm.df[,!colnames(comm.df) %in% c("id")]
  
  ##subset of the relatnedess rownames
  cophm<-coph[rownames(coph) %in% colnames(comm.df),]
  
  #trait distance matrix
  traitm<-traitdistance[rownames(traitdistance) %in% colnames(comm.df),]
  
#xytable
xytable<-xytab[id %in% rowsTocall]

  #Wrap the present
  #create unqiue filename for each rank
  fil<-paste("Data/",x,"Rank.RData",sep="")
  save(comm.df,Index_Space,cophm,traitm,xytable,file=fil)
  }
  }

finalize()
