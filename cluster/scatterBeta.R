#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals
#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW

##################


require(pbdMPI,quietly=TRUE)

init()

#Everyone say hello
comm.print(comm.rank(), all.rank = TRUE)
 
 #Everyone say hello
 
#Require libraries
  suppressMessages(require(reshape,quietly=TRUE,warn.conflicts=FALSE))
  suppressMessages(require(picante,quietly=TRUE,warn.conflicts=FALSE))
  suppressMessages(require(ggplot2,quietly=TRUE,warn.conflicts=FALSE))
  suppressMessages(require(parallel,quietly=TRUE,warn.conflicts=FALSE))
  suppressMessages(require(foreach,quietly=TRUE,warn.conflicts=FALSE))
  suppressMessages(require(GGally,quietly=TRUE,warn.conflicts=FALSE))
  suppressMessages(require(data.table,quietly=TRUE,warn.conflicts=FALSE))
  suppressMessages(require(vegan,quietly=TRUE,warn.conflicts=FALSE))
  
suppressMessages(source("Input/BrazilSourceFunctions.R"))

  ##If running locally set droppath
  
  droppath<-"/home1/02443/bw4sz/GlobalMammals/"
  
  setwd(droppath)


if(comm.rank()==0){

print("Rank 0 print")


 
   ###Define Source Functions, does this need to be run and distributed to all nodes, can they source simultaneously
  comm <- fread("Input/UniquesiteXspp.csv")    # as usual R read.table
  
  #read in xy data with original rows
  xytab<-fread("Output/xytable.csv")[,-1,with=F]
  setnames(xytab,colnames(xytab),c("V1","x","y","id"))
  
  print("siteXspp table:")
  print(comm[1:5,1:5,with=F])
  
  #make v1 a key column
  setkey(comm,id)
  
  print("Tables in memory:")
  tables()
  
  #Read in phylogeny
  tree<-read.tree("Input/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")
  
  #bring in traits
  traits <- read.table("Input/imputedmammals28apr14.txt",header=TRUE,row.names=1)
  
  #Read in cell by branch table made from source
  load(file="Input/tcellbr.RData")
  
  #Create all pairwise combinations of siteXspp
  z<-combn(comm$id,2)
  
  #index by key
  print(paste("Total number of iterations:",ncol(z)))
  
  #split rows into indices
  IndexFunction<-splitIndices(ncol(z),comm.size())
  print(paste("Length of IndexFunction is:",length(IndexFunction)))
  
  ###Send the subset matrix
  
  toScatter<-list()
for(x in 1:length(IndexFunction)){
    print(x)
    #Index of rows to call.
    Index_Space<-z[,IndexFunction[[x]]]
    
    #row positions
    rowsTocall<-unique(as.vector(Index_Space))
    comm.df<-data.frame(comm[id %in% rowsTocall])
    comm.df<-comm.df[,which(!apply(comm.df,2,sum)==0)]
    rownames(comm.df)<-comm.df$id
    comm.df<-comm.df[,!colnames(comm.df) %in% c("id")]
    
    ##subset of the tcellbr rownames 
    a<-tcellbr[rownames(tcellbr) %in% rowsTocall,]
    tcellbr_sub<-a[,which(!apply(a,2,sum)==0)]
    
rm(a)
print(lsos())
print(paste("Memory used:",mem()))

    #traits
    traitm<-traits[rownames(traits) %in% colnames(comm.df),]
    
    #xytable
    xytable<-xytab[id %in% rowsTocall]
    
    #Wrap the present
    #create unqiue filename for each rank

    toScatter[[x]]<-list(comm.df,Index_Space,tcellbr_sub,traitm,xytable)
  }

} else {
 toScatter<-NULL
xytab<-NULL
}

toScatterlist<-scatter(toScatter,rank.source=0)

print("scatter")

##################

###Define Source Functions

comm.df<-toScatterlist[[1]]
Index_Space<-toScatterlist[[2]]
tcellbr_sub<-toScatterlist[[3]]
traitm<-toScatterlist[[4]]
xytable<-toScatterlist[[5]]


timeF<-system.time(beta_out<-betaPar.scatter(comm.df,Index_Space,tcellbr_sub,traitm))

beta_out<-data.table(beta_out)

#Return timing argument to console
print(timeF)

#match to the xytable

# create all combinations of 2
# return as a data.table with these as columns `V1` and `V2`
# then count the numbers in each group
combs<-xytable[, rbindlist(combn(V1,2, 
                                 FUN = function(x) as.data.table(as.list(x)), simplify = F))]
setkey(combs,V1)
setkey(xytable,V1)

#combine both ways
combsV1<-xytable[combs]
setkey(combsV1,V2)
combsV2<-xytable[combsV1]

#set names
setnames(combsV2,colnames(combsV2),c("To.OriginalRow","To.xcord","To.ycord","To","From.OriginalRow","From.xcord","From.ycord","From"))

#Bring in pairwise betadiveristy data
print(head(beta_out))

#set keys to match
setkey(beta_out,To,From)

#set to character
combsV2[,To:=as.character(To)]
combsV2[,From:=as.character(From)]

setkey(combsV2,To,From)

#create combo columns
combR<-apply(beta_out,1,function(x){paste(sort(c(as.numeric(x[1]),as.numeric(x[2]))),collapse="_")})

combRxy<-apply(combsV2,1,function(x){paste(sort(c(as.numeric(x[4]),as.numeric(x[8]))),collapse="_")})

beta_out[,combo:=combR,]
combsV2[,combo:=combRxy]

setkey(beta_out,combo)
setkey(combsV2,combo)

head(mergeT<-merge(combsV2,beta_out))

print(dim(mergeT))

#try writing from all
comm.write.table(mergeT,"Output/FinalData.csv",row.names=FALSE)


finalize()

###############################
#Data Generation Complete
###############################