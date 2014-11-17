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
  print(A)
  print(B)
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
  
  #walk through each pair of cells and compute phylo betadiversity using MNNTDt
  MNdist <- sapply(rownames(comm),function(i){
    
    #set iterator
    A<-i
    
    #
    out<-lapply(rownames(comm)[1:(which(rownames(comm) == i))], function(B) {MNND(A,B,sp.list=sp.list,dists=dists)})
    names(out)<-rownames(comm)[1:(which(rownames(comm) == i))]
    return(out)
  })
  
  MNdist<-melt(MNdist)
  colnames(MNdist)<-c(nam,"To","From")
  return(MNdist)
  
}

##Compute Dimensions of betadiversity on all functions
beta_all<-function(comm,traitdist,coph){
  
  sp.list<-apply(comm,1,function(x){ 
    names(x[which(x==1)])
  })
  
  #taxonomic diversity
  tax<-taxF(comm)
  
  #name characters
  tax$To<-as.character(tax$To)
  tax$From<-as.character(tax$From)
  
  #phylogenetic
  phylo<-MNTDt(comm,coph,sp.list,nam="Phylo")
  
  #just use column names in the trait matrix
  comm.trait<-comm[,colnames(comm) %in% rownames(traitdist)]
  
  #trait species list
  sp.list.trait<-apply(comm.trait,1,function(x){ 
    names(x[which(x==1)])
  })
  
  #trait betadiversity
  trait<-MNTDt(comm=comm.trait,dists=traitdist,sp.list=sp.list.trait,nam="Trait")
  
  #merge together
  merge1<-merge(phylo,tax,by=c("To","From"))
  Allmetrics<-merge(merge1,trait,by=c("To","From"))
  
  #Combine with other metrics into one large dataframe
  return(Allmetrics)}

#Wrapper function!
betaPar.scatter<-function(toScatterMatrix,coph,traitdist){
    
  out<-beta_all(comm=toScatterMatrix,traitdist=traitdist,coph=coph)
  
  return(out)
}