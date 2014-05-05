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


#Branch matrix

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
matpsim <- function(com,tcellbr) # make sure nodes are labelled and that com and phyl species match
{
  
  # calculate full cell by cell phylobsim matrix
  
  psim <- lapply(rownames(tcellbr),function(j) {
    #print(j)
    cell_a <- j 
    
    out.cell<-lapply(rownames(tcellbr)[1:(which(rownames(tcellbr) == j))], 
                     nmatsim,cell_a=j)
    unlist(out.cell)
  })
  
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

nmatsim <- function(cell_b,cell_a,tcellbr) # samp = grid cell of interest
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
#Betadiversity functions
###############################

###############################
##Taxonomic Betadiversity

taxF<-function(comm){
  d<-as.matrix(vegdist(comm,binary=TRUE,upper=FALSE,diag=FALSE))
  sorenson<-melt(d)[melt(upper.tri(d))$value,]
  colnames(sorenson)<-c("To","From","Sorenson")
  return(sorenson)}

#################################
#Phylognetic Betadiversity
phyloF<-function(comm,tcellbr){
  #Compute cell matrix and melt it into a dataframe 
  betaSIM<-nmatsim(cell_a=rownames(comm)[1],cell_b=rownames(comm)[2],tcellbr)
  pmatSum<-data.frame(rownames(comm)[1],rownames(comm)[2],betaSIM)
  colnames(pmatSum)<-c("To","From","BetaSim")
  return(pmatSum)
}

#################################
#Trait Betadiversity

traitF<-function(comm,traits){
  
  #subset the data to find intersection of species
  #Trait frame needs to match siteXSpp table
  mon_cut<-traits[rownames(traits) %in% colnames(comm),]
  
  #species without traits, take them out for just this portion of the analysis, keep the assemblage lsit
  siteXspp_traits<-comm[,colnames(comm) %in% rownames(mon_cut)]
  
  #Data Checks
  if(!is.data.frame(siteXspp_traits)){
    melt.MNTD<-data.frame(MNTD=NA,To=rownames(comm)[1],From=rownames(comm)[2])
    return(melt.MNTD)
  }
     
     if( nrow(mon_cut)<2){
       melt.MNTD<-data.frame(MNTD=NA,To=rownames(comm)[1],From=rownames(comm)[2])
       return(melt.MNTD)
     }
     
     #if variance of trait column is 0, remove trait
     colvar<-!apply(mon_cut,2,var)==0
     mon_cut<-mon_cut[,colvar]
     prc_traits<-stats::prcomp(mon_cut,scale=TRUE)
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
     return(melt.MNTD)
}




##Compute Dimensions of betadiversity on all functions

#This is the workhorse function and could be split into parts
#The code computes phylogenetic (Phylosor), taxonomic (Sorenson) and trait betadiveristy (MNTD)
#for all sites entered into the list
#results are bound into a dataframe

#A phylogeny (tree), siteXspp matrix (comm) and trait matrix (traits) is required. 

beta_all<-function(comm=comm,tree=tree,traits=traits,tcellbr){
  
  #remove all species with no records in the row
  comm<-comm[,which(!apply(comm,2,sum)==0)]
  
  ##Taxonomic Betadiversity
  sorenson<-taxF(comm)
  
  ######Phylogenetic Betadiversity
  #Betasim from Holt 2013
  pmatSum<-phyloF(comm,tcellbr)
  
  #Merge with taxonomic
  Allmetrics0<-merge(pmatSum,sorenson,by=c("To","From"))
    
  melt.MNTD<-traitF(comm,traits)
  
  #Combine with other metrics into one large dataframe
  Allmetrics<-merge(Allmetrics0,melt.MNTD,by=c("To","From"))
  
  return(Allmetrics)}

#Wrapper function!


betaPar.scatter<-function(toScatterMatrix,toScatterIndex,tcellbr,rown,traits){
  
  ##############name the tcell matrix, the scatter function drops the names!!
  
  rownames(tcellbr)<-rown
  print(rownames(tcellbr)[1:10],all.rank=TRUE)
  
  #make sure the input data matches
  if(sum(!rownames(toScatterMatrix) %in% rownames(tcellbr))==0){
    print("Rownames of the matrix match the tcellbr")
  } else{
    
    stop("rownames of the matrix do not match the tcellbr")
    
  }
  
  #correct index
  if(sum(!rownames(toScatterMatrix) %in% unique(as.vector(toScatterIndex)))==0){
    print("Rownames of the matrix match the index of the matrix")
  } else{
    
    stop("rownames of the matrix do not match index of the matrix")
    
  }
  
  print(paste("Number of within loop calls:", ncol(toScatterIndex)))
  
  #Within a chunk, loop through the indexes and compute betadiversity
  holder<-apply(toScatterIndex,2,function(x) {
    print(x)
    #get the comm row
    comm.d<-toScatterMatrix[as.character(c(x[1],x[2])),]
    out<-beta_all(comm.d,tree=tree,traits=traits,tcellbr)
    return(out)
  }
  )
  
  #bind to a dataframe
  holder<-rbind.fill(holder)
  
  return(holder)
}
