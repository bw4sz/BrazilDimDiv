test1<-function(){
  
  #compute beta metrics
  #set flags
  phylosor.c<-FALSE
  beta.sim<-TRUE
  
  
  print(paste("Number of within loop calls:", ncol(toScatterIndex)))
  #Within a chunk, loop through the indexes and compute betadiversity
  holder<-apply(toScatterIndex,2,function(x) {
    print(x)
    #get the comm row
    comm.d<-toScatterMatrix[as.character(c(x[1],x[2])),]
    
    #if the rows are identical, betadiversity is 0
    if(duplicated(comm.d)[2]){
      out<-data.frame(To=rownames(comm.d)[[1]],From=rownames(comm.d)[[2]],BetaSim=0,Sorenson=0,MNTD=0)
    } else{
      
      out<-beta_all(comm.d,tree=tree,traits=traits,beta.sim,phylosor.c)
    }
    return(out)
  }
  )
  
  #bind to a dataframe
  holder<-rbind.fill(holder)
  
  return(holder)
}

Rprof(tmp <- tempfile())
testR<-test1()
Rprof()
summaryRprof(tmp)
unlink(tmp)