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
  compmat <- dists[J(Asp),Bsp,with=F]
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

MNTDt<-function(dists,sp.list,nam="value",ids){
  A<-as.character(ids[1])
  B<-as.character(ids[2])
  MNdist<-MNND(A,B,sp.list=sp.list,dists=dists)
  out<-data.frame(To=A,From=B,nam=MNdist)
  colnames(out)[3]<-nam
  return(out)
}

##Compute Dimensions of betadiversity on all functions
beta_all<-function(comm,traitdist,coph,sp.list.phylo,sp.list.trait,ids){
  
  #taxonomic diversity
  cm<-melt(sp.list.phylo[as.character(ids)])
  cm$PA<-1
  tax<-taxF(dcast(cm,L1~value,value.var='PA',fill = 0)[,-1])[2,]
  
  #name characters
  tax$To<-as.character(ids[1])
  tax$From<-as.character(ids[2])
  
  #phylogenetic
  phylo<-MNTDt(dists=coph,sp.list=sp.list.phylo,nam="Phylo",ids)
  
  #trait betadiversity
  trait<-MNTDt(dists=traitdist,sp.list=sp.list.trait,nam="Trait",ids)
  
  #merge together
  merge1<-merge(phylo,tax,by=c("To","From"))
  Allmetrics<-merge(merge1,trait,by=c("To","From"))
  
  #Combine with other metrics into one large dataframe
  return(Allmetrics)}

#Wrapper function!
betaPar.scatter<-function(toScatterIndex,coph,traitdist,sp.list.phylo,sp.list.trait){
    
  #Within a chunk, loop through the indexes and compute betadiversity
  holder<-apply(toScatterIndex,2,function(x) {
    #get the comm row
    out<-beta_all(coph=coph,traitdist=traitdist,sp.list.phylo = sp.list.phylo,sp.list.trait = sp.list.trait,ids=x)
    return(out)
  })
  
  #bind to a dataframe
  holder<-rbind.fill(holder)  
  return(holder)
}

#Function to randomly sample species in each assemblage
sampleS<-function(x,ident){
  row<-ident[x,]
  #Get species for assemblage A
  A<-sample(splist,size=row[[1]],replace=F)
  # Get species for Assemblage B
  B<-sample(splist,size=row[[2]],replace=F)
  assemb<-list(A,B)
  names(assemb)<-c("A","B")

  #Format the data frame
  assemblages<-melt(assemb)
  assemblages<-as.data.frame.array(t(table(assemblages)))
  return(assemblages)
}