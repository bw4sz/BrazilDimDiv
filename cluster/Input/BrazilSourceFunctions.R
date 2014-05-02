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
#Global Betadiversity function
###############################

#This is the workhorse function and could be split into parts
#The code computes phylogenetic (Phylosor), taxonomic (Sorenson) and trait betadiveristy (MNTD)
#for all sites entered into the list
#results are bound into a dataframe

#A phylogeny (tree), siteXspp matrix (comm) and trait matrix (traits) is required. 

beta_all<-function(comm=comm,tree=tree,traits=traits,tcellbr){
  
  
  if(sum(comm)==0){return(NA)}
  
  #remove all species with no records in the row
  comm<-comm[,which(!apply(comm,2,sum)==0)]
  
  #####################################
  ##Taxonomic Betadiversity
  
  d<-as.matrix(vegdist(comm,binary=TRUE,upper=FALSE,diag=FALSE))
  sorenson<-melt(d)[melt(upper.tri(d))$value,]
  colnames(sorenson)<-c("To","From","Sorenson")
  ######################################
   
  ######Phylogenetic Betadiversity
  #Betasim from Holt 2013
  
  #I broke this into seperate functions since it doesn't need to happen on each cell
  
  #From the initial rank 0 i computed the branching matrix and stored it as branch.out
  
    #Compute cell matrix and melt it into a dataframe 
    betaSIM<-nmatsim(cell_a=rownames(comm)[1],cell_b=rownames(comm)[2],tcellbr)
    pmatSum<-data.frame(rownames(comm)[1],rownames(comm)[2],betaSIM)
    colnames(pmatSum)<-c("To","From","BetaSim")
    
    #Merge with taxonomic
    Allmetrics0<-merge(pmatSum,sorenson,by=c("To","From"))
  
  
  #Trait frame needs to match siteXSpp table
  mon_cut<-traits[rownames(traits) %in% colnames(comm),]
  
  #species without traits, take them out for just this portion of the analysis, keep the assemblage lsit
  siteXspp_traits<-comm[,colnames(comm) %in% rownames(mon_cut)]
  
  #direct traits, need to be standardized
 
  
  #If you wanted to collapse traits into PC axis using: 

#if variance is 0, remove trait
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
  
  #Combine with other metrics into one large dataframe
  Allmetrics<-merge(Allmetrics0,melt.MNTD,by=c("To","From"))
  
  return(Allmetrics)}


#Wrapper function!


betaPar.scatter<-function(toScatterMatrix,toScatterIndex,tcellbr,rown){
  
##############name the tcell matrix, the scatter function drops the names!!

rownames(tcellbr)<-rown
comm.print(rownames(tcellbr)[1:10],all.rank=TRUE)


#make sure there is matches
if(sum(!rownames(toScatterMatrix) %in% rownames(tcellbr))==0){
print("Rownames of the matrix match the tcellbr")
} else{
print("rownames of the matrix DO not match the tcellbr")
break
}

#correct index
if(sum(!rownames(toScatterMatrix) %in% unique(as.vector(toScatterIndex)))==0){
print("Rownames of the matrix match the index of the matrix")
} else{
print("rownames of the matrix DO not match index of the matrix")
break
}
    
  print(paste("Number of within loop calls:", ncol(toScatterIndex)))

  #Within a chunk, loop through the indexes and compute betadiversity
  holder<-apply(toScatterIndex,2,function(x) {

    #get the comm row
    comm.d<-toScatterMatrix[as.character(c(x[1],x[2])),]
    
    #if the rows are identical, betadiversity is 0
    if(duplicated(comm.d)[2]){
      out<-data.frame(To=rownames(comm.d)[1],From=rownames(comm.d)[2],BetaSim=0,Sorenson=0,MNTD=0)
    } else{
      
      out<-beta_all(comm.d,tree=tree,traits=traits,tcellbr)
    }
    return(out)
  }
  )
  
  #bind to a dataframe
  holder<-rbind.fill(holder)
  
  return(holder)
}
