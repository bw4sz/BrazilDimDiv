#Taxonomic, Phylogenetic and Functional Diversity in South American Hummingbirds
#Manuscript submitted to the American Naturalist

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW
#Boris A. Tinoco, Juan L. Parra, PhD, Leone M. Brown, PhD, Gary Stiles, PhD, Jim A. McGuire, PhD, Catherine H. Graham, PhD


#############
#Set testing flag, 
#this subsitutes the real data for a tiny dummy example for debugging

testing<-FALSE

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

#Everyone say hello
#comm.print(comm.rank(), all.rank = TRUE)

##If running locally set droppath

droppath<-"/home1/02443/bw4sz/GlobalMammals/"

setwd(droppath)

###Define Source Functions, does this need to be run and distributed to all nodes, can they source simultaneously

suppressMessages(source("Input/BrazilSourceFunctions.R"))


########################################
###############Read in data on the node
########################################

#read in all same ranks
#siteXspp <- comm.read.csv("Input/siteXspp1dgr.csv",read.method="common")

##barrier method
print("barrier method start")
for(i.rank in 0:(comm.size() - 1)){
  if(i.rank == comm.rank()){
    	
	siteXspp <- read.csv("Input/siteXspp1dgr.csv",row.names=1)    # as usual R read.table
	
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
  }
  barrier()
}

print("data loaded!")

  #Remove xy data
  siteXspp<-siteXspp[,!colnames(siteXspp) %in% c("X","x","y","rich")]
  
#Remove lines with less than 2 species
  richness<-apply(siteXspp,1,sum)
  keep<-which(richness > 1)
  siteXspp<-siteXspp[keep,]
  #Get entire species list
  splist<-colnames(siteXspp)

#remove species in the siteXspp that are not in phylogeny
siteXspp<-siteXspp[,colnames(siteXspp) %in% tree$tip.label]
print(dim(siteXspp))

#Do we want to subset the data for a test case? 
#If so, add which rows below
comm<-siteXspp[1:1000,]

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

