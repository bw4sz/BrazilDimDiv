#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals
#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW


##################

library(pbdMPI,quietly=TRUE,warn.conflicts=FALSE,verbose=F)

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
comm <- fread("Input/UniquesiteXspp.csv")    # as usual R read.table
#make v1 a key column
setkey(comm,id)

#read in xy data with original rows
xytab<-fread("Output/xytable.csv")[,-1,with=F]
setnames(xytab,colnames(xytab),c("x","y","V1","id"))

print("siteXspp table:")
print(comm[1:5,1:5,with=F])


print("Tables in memory:")
tables()

#Bring in relatedness matrix
coph<-fread("Input/cophenetic.csv")
setkey(coph,V1)
setcolorder(coph,c("V1",coph$V1))

#bring in trait distance matrix
traitdistance<-fread("Input/traitdistanceLog.csv")
setkey(traitdistance,V1)

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
  comm.df<-comm[J(rowsTocall)]
  setkey(comm.df,id)
  
  #remove species
  c2<-names(which(comm.df[,lapply(.SD,sum)==0]))
  
  #subset columns
  comm.df<-comm.df[,!c(c2),with=F]
    
  ##subset of the relatnedess rownames
  coph_keep<-colnames(coph)[colnames(coph) %in% colnames(comm.df)]
  cophm<-coph[J(colnames(comm.df)),c("V1",coph_keep),with=F]
  cophm<-cophm[-1,]
  setkey(cophm,V1)
  
  #trait distance matrix
  trait_keep<-colnames(traitdistance)[colnames(traitdistance) %in% colnames(comm.df)]
  traitm<-traitdistance[J(trait_keep),c("V1",trait_keep),with=F]
  setkey(traitm,V1)

  #xytable
xytable<-xytab[id %in% rowsTocall]

  #Wrap the present
  #create unqiue filename for each rank
  fil<-paste("Data/",x,"Rank.RData",sep="")
  save(comm.df,Index_Space,cophm,traitm,xytable,file=fil)
  }
  }

finalize()
