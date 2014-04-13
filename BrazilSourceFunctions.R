#Source functions for Global mammal diversity metrics

##Append a position to a list
lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}


##Get species list from cell number
#This helpful function get species list, but make sure commid is in ""

Getsplist<-function(commID){
  names(siteXspp[commID,which(siteXspp[commID,]==1)])
}


############Trait Betadiversity Function,
#Computes the mean nearest neighbor

MNND <- function(A,B,sp.list,dists)
{
  
  Asp     <- sp.list[[A]]
  Bsp     <- sp.list[[B]]
  compmat <- dists[Asp,Bsp]
  Ann     <- apply(as.matrix(compmat),1,min)
  Bnn     <- apply(as.matrix(compmat),2,min)
  Dnn     <- mean(c(Ann, Bnn))
  #turn    <- min(c(mean(Ann),mean(Bnn)))
  #nest    <- Dnn - turn
  #res <- c(Dnn,turn,nest)
  res<-c(Dnn)
  #names(res) <- c("MNND","MNNDturn","MNNDnest")
  names(res) <- c("MNND")
  
  res
}



##Ben Holt's matrix phylo betadiversity function, let's break this into pieces, so its not-redundant on each call

#Function gets the edge matrix for all nodes, takes awhile
branching<-function(phyl){
  require(phylobase)  # detail information for all phylo branches
  new <- phylo4(phyl)
  dat <- as.data.frame(print(new))
  allbr <- dat$edge.length
  names(allbr) <- getEdge(new)
  return(list(dat,new,allbr))
}


#Function computes the cell by branch matrix

psimbranches<-function(phyl=tree,com=comm,branch.out=branch.out){
  #unpack the output from the branching function
  dat<-branch.out[[1]]
  new<-branch.out[[2]]
  allbr<-branch.out[[3]]
  
  #get species list
  spp <- colnames(com)
  
  # create a list of phy branches for each species
  brs <-  foreach(i = spp, .packages = "phylobase") %do% #this loop makes a list of branches for each species
{  
  print(which(spp == i)/length(spp))
  print(date())
  print(i)
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
  
  #print("brs")
  
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
  
  #print("spp_br")
  
  spp_br <<- spp_br
  
  tcellbr <- foreach(j = rownames(com), .combine = "rbind") %do% {
    cellbr(j,spp_br,com)
  }
  
  #name columns
  rownames(tcellbr) <- rownames(com)
  tcellbr <<- tcellbr
  
  #print("cell_br")
  #print(tcellbr)
  return(tcellbr)
}

#create a species by phy branch matrix, broken out from original function

matpsim <- function(tcellbr) # make sure nodes are labelled and that com and phyl species match
{
  # calculate full cell by cell phylobsim matrix
  
  psim <- foreach(j = rownames(tcellbr)) %do% {
    #print(j)
    cell_a <- j 
    unlist(
      lapply(rownames(tcellbr)[1:(which(rownames(tcellbr) == j))], 
             nmatsim,cell_a=j
      )
    )
  }
  
  #print("psim")
  psim <- do.call("rbind", psim)
  rownames(psim) <- rownames(tcellbr)
  colnames(psim) <- rownames(tcellbr)    
  psim[upper.tri(psim)] <- psim[lower.tri(psim)]
  #psim <- as.dist(psim)
  #print("its over")
  return(psim)
}

#calculate phylosim between cells, broken out from original function

nmatsim <- function(cell_b,cell_a) # samp = grid cell of interest
{
  a_br  <- tcellbr[cell_a,]
  b_br <- tcellbr[cell_b,]
  s_br <- rbind(a_br,b_br)
  s_br <- as.matrix(s_br[,colSums(s_br) >0])
  pa_br <- s_br > 0
  both <- s_br[1,colSums(pa_br > 0)==2]
  ubr <- s_br[,colSums(pa_br > 0)==1]
  a_ubr <- as.matrix(ubr)[1,]
  b_ubr <- as.matrix(ubr)[2,]
  psim <- 1 - (sum(both,na.rm=T)/(min(sum(a_ubr,na.rm=T),sum(b_ubr,na.rm=T))+sum(both,na.rm=T)))
  return(psim)
}


# function to give the pres/abs of phy branches withi cell i, broken out from original function
cellbr <- function(i,spp_br, com)
{
  i_spp <- rownames(spp_br)[com[i,]>0]
  if(length(i_spp) > 1)
  {i_br  <- as.numeric(apply(spp_br[i_spp,],2,max))}
  else  {i_br  <- as.numeric(spp_br[i_spp,]) }
  names(i_br) <- colnames(spp_br)
  return(i_br)
}

###############################
#Global Betadiversity function
###############################

#This is the workhorse function and could be split into parts
#The code computes phylogenetic (Phylosor), taxonomic (Sorenson) and trait betadiveristy (MNTD)
#for all sites entered into the list
#results are bound into a dataframe

#A phylogeny (tree), siteXspp matrix (comm) and trait matrix (traits) is required. 

beta_all<-function(comm=comm,tree=tree,traits=traits,tcellbr=tcellbr,beta.sim,phylosor.c){
  
  if(sum(comm)==0){return(NA)}
  
  #remove all species with no records in the row
  comm<-comm[,which(!apply(comm,2,sum)==0)]
  
  #####################################
  ##Taxonomic Betadiversity
  
  d<-as.matrix(vegdist(comm,binary=TRUE,upper=FALSE,diag=FALSE))
  sorenson<-melt(d)[melt(upper.tri(d))$value,]
  colnames(sorenson)<-c("To","From","Sorenson")
  ######################################
  
  if(phylosor.c==TRUE){
    ######################################
    #Phylogenetic Betadiversity
    
    #Phylosor Calculation see Bryant 2008
    phylo.matrix<-as.matrix(phylosor(comm,tree))
    diag(phylo.matrix)<-NA
    Phylosor.phylo<-melt(phylo.matrix)
    colnames(Phylosor.phylo)<-c("To","From","Phylosor.Phylo")
    Phylosor.phylo$Phylosor.Phylo<-1-Phylosor.phylo$Phylosor.Phylo
    Allmetrics0<-merge(Phylosor.phylo,sorenson,by=c("To","From"))
  }
  #######################################
  
  ######Phylogenetic Betadiversity
  #Betasim from Holt 2013
  
  ###Phylosor is really slow from the PD function, try ben holt's betasim
  #I broke this into seperate functions since it doesn't need to happen on each cell
  
  #From the initial rank 0 i computed the branching matrix and stored it as branch.out
  if(beta.sim==TRUE){
    
    #Compute cell matrix and melt it into a dataframe 
    betaSIM<-matpsim(tcellbr)
    pmatSum<-melt(betaSIM)
    colnames(pmatSum)<-c("To","From","BetaSim")
    
    #Merge with taxonomic
    Allmetrics0<-merge(pmatSum,sorenson,by=c("To","From"))
  }
  
  #Trait frame needs to match siteXSpp table
  mon_cut<-traits[rownames(traits) %in% colnames(comm),]
  
  #species without traits, take them out for just this portion of the analysis, keep the assemblage lsit
  siteXspp_traits<-comm[,colnames(comm) %in% rownames(mon_cut)]
  
  #direct traits, need to be standardized
  
  #   #Zscores, standardized by sd and subtracted means
  #   means<-apply(mon_cut,2,mean)
  #   Bill<-mon_cut$Bill - means["Bill"]/sd(mon_cut$Bill)
  #   Mass<-mon_cut$Mass - means["Mass"]/sd(mon_cut$Mass)
  #   WingChord<-(mon_cut$WingChord - means["WingChord"])/sd(mon_cut$WingChord)
  #   z.scores<-data.frame(Bill,Mass,WingChord)
  #   rownames(z.scores)<-rownames(mon_cut)
  #   newSGdist<-dist(z.scores,method="euclidean")
  #   
  
  #If you wanted to collapse traits into PC axis using: 
  prc_traits<-prcomp(mon_cut,scale=TRUE)
  newSGdist <- dist(prc_traits$x)
  
  #create sp.list for trait function
  sp.list<-lapply(rownames(siteXspp_traits),function(k){
    x<-siteXspp_traits[k,]
    names(x[which(x==1)])
  })
  
  names(sp.list)<-rownames(siteXspp_traits)
  dists <- as.matrix(newSGdist)
  
  rownames(dists) <- rownames(mon_cut)
  colnames(dists) <- rownames(mon_cut)
  
  #walk through each pair of cells and compute trait betadiversity using MNNTD
  sgtraitMNTD <- sapply(rownames(siteXspp_traits),function(i){
    
    #Iterator count
    #print(round(which(rownames(siteXspp_traits)==i)/nrow(siteXspp_traits),3))
    
    #set iterator
    A<-i
    
    #
    out<-lapply(rownames(siteXspp_traits)[1:(which(rownames(siteXspp_traits) == i))], function(B) {MNND(A,B,sp.list=sp.list,dists=dists)})
    names(out)<-rownames(siteXspp_traits)[1:(which(rownames(siteXspp_traits) == i))]
    return(out)
  })
  
  #name trait matrix
  names(sgtraitMNTD) <- rownames(siteXspp_traits)
  melt.MNTD<-melt(sgtraitMNTD)
  colnames(melt.MNTD)<-c("MNTD","To","From")
  
  #Combine with other metrics into one large dataframe
  Allmetrics<-merge(Allmetrics0,melt.MNTD,by=c("To","From"))
  
  return(Allmetrics)}


#Wrapper function!

#betaPar function takes in the siteXspp argument (comm), the rank of the node (rank), and 
#the chunks for which to split indices. The function finds all pairwise comparisons of an index (slow)
#breaks the index into chunk pieces, and runs these chunks on seperate nodes of the cluster

betaPar<-function(comm,rankNumber,chunks,beta.sim,phylosor.c){
  #Create all pairwise combinations of siteXspp
  z<-combn(nrow(comm),2)
  
  #split rows into indices, we want each loop to take about an hour, 
  #THe function is initially timed at 20 seconds, 
  IndexFunction<-splitIndices(ncol(z),chunks)
  
  ###Divide the indexes, ########THE ONE IS CRUCIAL HERE< THIS NEEDS TO BE RANKED ON PBDMPI
  Index_Space<-z[,IndexFunction[[rankNumber]]]
  
  #Create an output to hold function container, there are as many rows as columns in the combinations, and there are five columns for the output data
  holder<-list()
  
  #Within a chunk, loop through the indexes and compute betadiversity
  for (x in 1:ncol(Index_Space)){
    print(x)
    #Grab the correct index
    index_col<-Index_Space[,x] 
    
    #get the comm row
    comm.d<-comm[c(index_col[[1]],index_col[[2]]),]
    
    #if the rows are identical, betadiversity is 0
    if(sum(!comm.d[1,]==comm.d[2,])==0){
      holder[[x]]<-data.frame(To=index_col[[1]],From=index_col[[2]],BetaSim=0,Sorenson=0,MNTD=0)
    } else{
      
      #compute beta metrics
      out<-beta_all(comm.d,tree=tree,traits=traits,tcellbr,phylosor.c,beta.sim)
      holder[[x]]<-out
    }
  }
  
  #bind to a dataframe
  holder<-rbind.fill(holder)
  
  return(holder)}
