#Source functions for DimDiv

lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

#define a helpful function get species list commid MUST be in quotes!

Getsplist<-function(commID){
  names(siteXspp[commID,which(siteXspp[commID,]==1)])
}

#Define Betadiversity function

beta_all<-function(comm,tree,traits){
  
  #####################################
  ##Taxonomic Betadiversity
  d<-as.matrix(vegdist(comm,binary=TRUE,upper=FALSE,diag=FALSE))
  sorenson<-melt(d)[melt(upper.tri(d))$value,]
  colnames(sorenson)<-c("To","From","Sorenson")
  ######################################
  
  #Phylosor Calculation see Bryant 2008
  phylo.matrix<-as.matrix(phylosor(comm,tree))
  diag(phylo.matrix)<-NA
  Phylosor.phylo<-melt(phylo.matrix)
  colnames(Phylosor.phylo)<-c("To","From","Phylosor.Phylo")
  Phylosor.phylo$Phylosor.Phylo<-1-Phylosor.phylo$Phylosor.Phylo
  
  #Trait frame needs to match siteXSpp table
  mon_cut<-traits[rownames(traits) %in% colnames(comm),]
  
  #There are eight species without traits, take them out for just this portion of the analysis, keep the assemblage lsit
  siteXspp_traits<-comm[,colnames(comm) %in% rownames(mon_cut)]
  
  #   #Zscores, standardized by sd and subtracted means
  #   means<-apply(mon_cut,2,mean)
  #   Bill<-mon_cut$Bill - means["Bill"]/sd(mon_cut$Bill)
  #   Mass<-mon_cut$Mass - means["Mass"]/sd(mon_cut$Mass)
  #   WingChord<-(mon_cut$WingChord - means["WingChord"])/sd(mon_cut$WingChord)
  #   z.scores<-data.frame(Bill,Mass,WingChord)
  #   rownames(z.scores)<-rownames(mon_cut)
  #   newSGdist<-dist(z.scores,method="euclidean")
  #   
  #If you wanted a PCA 
  prc_traits<-prcomp(mon_cut,scale=TRUE)
  newSGdist <- dist(prc_traits$x)
  
  #create sp.list
  sp.list<-lapply(rownames(siteXspp_traits),function(k){
    x<-siteXspp_traits[k,]
    names(x[which(x==1)])
  })
  
  names(sp.list)<-rownames(siteXspp_traits)
  dists <- as.matrix(newSGdist)
  
  rownames(dists) <- rownames(mon_cut)
  colnames(dists) <- rownames(mon_cut)
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
  
  names(sgtraitMNTD) <- rownames(siteXspp_traits)
  melt.MNTD<-melt(sgtraitMNTD)
  colnames(melt.MNTD)<-c("MNTD","To","From")
  
  #Combine with other metrics
  Allmetrics0<-merge(Phylosor.phylo,melt.MNTD,by=c("To","From"))
  Allmetrics<-merge(Allmetrics0,sorenson,by=c("To","From"))
  
  return(Allmetrics)}


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

#Function to take in a apir of

betaPar<-function(comm){
  
#Create all pairwise combinations of siteXspp
z<-expand.grid(1:nrow(comm),1:nrow(comm))

#make a unique key
z$key <- apply(z, 1, function(x)paste(sort(x), collapse=''))

#get rid of duplicates
expandd<-subset(z, !duplicated(z$key))

colnames(expandd)<-c("To","From","Key")

#get rid of diagonal
expandd<-expandd[!expandd$To==expandd$From,1:2]

#split rows into indices.
IndexFunction<-splitIndices(nrow(expandd),2)

###Divide the indexes, ########THE ONE IS CRUCIAL HERE< THIS NEEDS TO BE RANKED ON PBDMPI
Index_Space<-expandd[IndexFunction[[1]],]

#for each combination

for (x in 1:nrow(IndexSpace)){
 index_row<-Index_Space[x,] 

 #get the comm row
 comm.d<-comm[c(index_row$To,index_row$From),]
 
 out<-beta_all(comm.d,tree=tree,traits=traits)
 return(out)
}
