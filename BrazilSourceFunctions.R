#Source functions for DimDiv


lappend <- function(lst, obj) {
  lst[[length(lst)+1]] <- obj
  return(lst)
}

#define a helpful function get species list commid MUST be in quotes!

Getsplist<-function(commID){
  names(siteXspp[commID,which(siteXspp[commID,]==1)])
}
