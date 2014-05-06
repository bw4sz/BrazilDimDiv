#Compute Phylogenetic Cell by Branch matrix for all unique comparisons using Betasim

suppressMessages(require(picante,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(foreach,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(doSNOW,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(phylobase,quietly=TRUE,warn.conflicts=FALSE))

#Set dropbox path

droppath<-"/home1/02443/bw4sz/GlobalMammals/"

setwd(droppath)

####Read in data
#Read in species matrix

siteXspp <- as.matrix(read.csv("Output/UniquesiteXspp.csv"))

print(siteXspp[1:5,1:5])

rownames(siteXspp)<-siteXspp[,"id"]

#Just get the species data
siteXspp<-siteXspp[,!colnames(siteXspp) %in% c("X","x","y","rich","V1","V0","id")]

#read in phylogeny

tree<-read.tree("Input/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")

print(tree)

matpsim <- function(phyl, com) # make sure nodes are labelled and that com and phyl species match
{
  
  require(phylobase)  # detail information for all phylo branches
  new <- phylo4(phyl)
  dat <- as.data.frame(print(new))
  allbr <- dat$edge.length
  names(allbr) <- getEdge(new)
  spp <- colnames(com)

  
  cl <- getMPIcluster() # create parellel clusters
  registerDoSNOW(cl)
  
  # create a list of phy branches for each species
  
  brs <-  foreach(i = spp, .packages = "phylobase") %dopar% #this loop makes a list of branches for each species
{  

  brsp <- vector()
  br   <- as.numeric(rownames(dat[which(dat$label==i),]))
  repeat{
    brsn <- getEdge(new,br)
    brsl <- dat[names(brsn),"edge.length"]
    names(brsl) <- brsn
    brsp <- c(brsp, brsl)
    br   <- dat[br,3]
    if(br == 0) {break}
  }
  
  brsp
}
  names(brs) <- spp
  
    
  print("brs")
  
  # create a species by phy branch matrix
  
  spp_br <- matrix(0,nrow = length(spp), ncol = length(allbr))
  rownames(spp_br) <- spp
  colnames(spp_br) <- names(allbr)
  
  for(i in spp)
  {
    spp_br[i,names(brs[[i]])] <- brs[[i]]
  }
  
  spp_br <- spp_br[,-(ncol(com)+1)] # removes root
  spp_br <- spp_br[,!colSums(spp_br) %in% c(0,nrow(spp_br))]  # here take out all common branches instead
  
  print("spp_br")
  
  spp_br <<- spp_br
  
  # function to give the pres/abs of phy branches withi cell i
  
  cellbr <- function(i,spp_br, com)
  {
    i_spp <- rownames(spp_br)[com[i,]>0]
    if(length(i_spp) > 1)
    {i_br  <- as.numeric(apply(spp_br[i_spp,],2,max))}
    else  {i_br  <- as.numeric(spp_br[i_spp,]) }
    names(i_br) <- colnames(spp_br)
    return(i_br)
  }
  
  print(rownames(com)[1:10])
  
  tcellbr <- foreach(j = rownames(com), .combine = "rbind") %dopar% {cellbr(j,spp_br,com)}
  
  print("cell_br")
  rownames(tcellbr) <- rownames(com)
  tcellbr <<- tcellbr
return(tcellbr)
}

#compute
tcellbr<-matpsim(phyl=tree,com=siteXspp)

print("function computed")

print(tcellbr[1:5,1:5])
print(dim(tcellbr))


#Write to file
save(tcellbr,file="Output/tcellbr.RData")

print("imaged")
