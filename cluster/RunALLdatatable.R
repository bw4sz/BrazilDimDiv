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

droppath<-"C:/Users/Jorge/Documents/BrazilDimDiv/cluster"

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
system.time(comm<-siteXspp[V1 <= 100])

print(paste("There are in comm NAs:",sum(is.na(comm))))

#Create all pairwise combinations of siteXspp
z<-combn(comm$V1,2)

#index by key

#print(paste("Total number of iterations:",ncol(z)))

#split rows into indices, we want each loop to take about an hour, 
#THe function is initially timed at 20 seconds, 
IndexFunction<-splitIndices(ncol(z),chunks)

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
}

#####################################################
##############Compute Betadiversity##################
#####################################################

comm.print(paste("Total nodes is:", comm.size()))

.rank<-comm.rank() + 1
timeF<-system.time(beta_out<-betaPar(comm,rankNumber=.rank,chunks=comm.size(),phylosor.c=FALSE,beta.sim=TRUE))

#combine rownames with xycoordinates


#Test alterantive methods
#benchmark the two approaches across 10 runs.
#require(rbenchmark)
#a<-benchmark(replications=1,
#betaPar(comm,1,2,beta.sim=TRUE,phylosor.c=FALSE),
#betaPar(comm,1,2,beta.sim=FALSE,phylosor.c=TRUE))

#Correlation among dimensions
#ggpairs(beta_out[,-c(1,2)],cor=TRUE)

#Return timing argument to console
print(timeF)

#try writing from all
comm.write.table(beta_out,"beta_out.csv")

finalize()

#Data Generation Complete
##########################################################################################
##########################################################################################
