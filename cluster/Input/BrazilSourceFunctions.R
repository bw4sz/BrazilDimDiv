#Source functions for Global mammal diversity metrics

#Memory function
# improved list of objects

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
  turn    <- min(c(mean(Ann),mean(Bnn)))
  #nest    <- Dnn - turn
  #res <- c(Dnn,turn,nest)
  res<-c(turn)
  #names(res) <- c("MNND","MNNDturn","MNNDnest")
  names(res) <- c("MNNDt")
  
  res
}


###############################
#Betadiversity functions
###############################

###############################
##Taxonomic Betadiversity

taxF<-function(comm){
  d<-betapart.core(comm)
  bsim<-d$min.not.shared/(d$min.not.shared + d$shared)
  rownames(bsim)<-rownames(comm)
  colnames(bsim)<-rownames(comm)
  bsim<-melt(bsim)
  colnames(bsim)<-c("To","From","Tax")
  return(bsim)
}


#################################
#Mean Taxon Betadiversity, turnover only
####################################

MNTDt<-function(comm,dists,sp.list,nam="value"){
  A<-rownames(comm)[1]
  B<-rownames(comm)[2]
  MNdist<-MNND(A,B,sp.list=sp.list,dists=dists)
  out<-data.frame(To=A,From=B,nam=MNdist)
  colnames(out)[3]<-nam
  return(out)
}

##Compute Dimensions of betadiversity on all functions
beta_all<-function(comm,traitdist,coph,sp.list.phylo,sp.list.trait){
  
  #taxonomic diversity
  tax<-taxF(comm)
  
  #name characters
  tax$To<-as.character(tax$To)
  tax$From<-as.character(tax$From)
  
  #phylogenetic
  phylo<-MNTDt(comm,coph,sp.list.phylo,nam="Phylo")
  
  #just use column names in the trait matrix
  comm.trait<-comm[,colnames(comm) %in% rownames(traitdist)]
  
  #trait betadiversity
  trait<-MNTDt(comm=comm.trait,dists=traitdist,sp.list=sp.list.trait,nam="Trait")
  
  #merge together
  merge1<-merge(phylo,tax,by=c("To","From"))
  Allmetrics<-merge(merge1,trait,by=c("To","From"))
  
  #Combine with other metrics into one large dataframe
  return(Allmetrics)}

#Wrapper function!
betaPar.scatter<-function(toScatterMatrix,toScatterIndex,coph,traitdist,sp.list.phylo,sp.list.trait){
    
  #Within a chunk, loop through the indexes and compute betadiversity
  holder<-apply(toScatterIndex,2,function(x) {
    #get the comm row
    comm.d<-toScatterMatrix[as.character(c(x[1],x[2])),]
    out<-beta_all(comm = comm.d,coph=cophm,traitdist=traitm,sp.list.phylo = sp.list.phylo,sp.list.trait = sp.list.trait)
    return(out)
  })
  
  #bind to a dataframe
  holder<-rbind.fill(holder)  
  return(holder)
}