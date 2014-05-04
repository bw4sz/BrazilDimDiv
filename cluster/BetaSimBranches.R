#Define Function



suppressMessages(require(picante,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(foreach,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(doSNOW,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(phylobase,quietly=TRUE,warn.conflicts=FALSE))

#Set dropbox path

droppath<-"/home1/02443/bw4sz/GlobalMammals/"

setwd(droppath)

####Read in data
#Read in species matrix

siteXspp <- read.csv("Input/UniquesiteXspp.csv")

print(siteXspp[1:5,1:5])

#Just get the species data, starts on column 33 for this example
siteXspp<-siteXspp[,!colnames(siteXspp) %in% c("X","x","y","rich","V1","V0")]

#Remove lines with less than 2 species
richness<-apply(siteXspp,1,sum)
keep<-which(richness > 1)
siteXspp<-siteXspp[keep,]

#Get entire species list
splist<-colnames(siteXspp)

#Read in phylogeny
tree<-read.tree("Input/Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")

#remove species in the siteXspp that are not in phylogeny
siteXspp<-siteXspp[,colnames(siteXspp) %in% tree$tip.label]

print(dim(siteXspp))

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
  print(which(spp == i)/length(spp))
  print(date())
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
  
    
  tcellbr <- foreach(j = rownames(com), .combine = "rbind") %dopar% {cellbr(j,spp_br,com)}
  
  print("cell_br")
  rownames(tcellbr) <- rownames(com)
  tcellbr <<- tcellbr

}

#compute
tcellbr<-matpsim(phyl=tree,com=siteXspp)

print("function computed")
print(dim(tcellbr))


#Write to file
save(tcellbr,file="Input/tcellbr.RData")

print("imaged")

