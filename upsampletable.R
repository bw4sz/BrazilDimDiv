###Create a siteXspp table at a new resolution 

require(reshape)
require(foreach)
require(doSNOW)
require(raster)
require(data.table)

##read in data
droppath<-"C:/Users/Jorge/Dropbox"

siteXspp<-fread(paste(droppath,"Dimensions Data/siteXspp.csv",sep="/"),nrow=35000)

#set the data table as a key
setnames(siteXspp,"V1","XY")
setkey(siteXspp,XY)

#create a molleweide raster with the desired resolution
a<-raster()
res(a)<-1

#project raster
a.proj<-projectRaster(a,crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")

#Get XY coordinates for the XY position
XY<-colsplit(siteXspp$XY," ",c("X","Y"))

#get the cell number from the projected raster
cellnumber<-cellFromXY(a.proj,XY)

#create an index of XY to cell number in old and new table
sp.index<-data.frame(XY=siteXspp$XY,cellnumber)

#make into a new data.table
sp.indexdt<-data.table(sp.index)
setkey(sp.indexdt,"XY")

#bind the new cellnumbers to the old matrix
siteXsppM<-merge(sp.indexdt,siteXspp)

tables()

foo<-function(x){
  (sum(x) >= 1)*1
} 

#aggregate by cellnumber on the new merged data.table
siteXsppSUM<-siteXsppM[,lapply(.SD,foo),by=cellnumber,.SDcols = colnames(siteXsppM)[!colnames(siteXsppM) %in% c("cellnumber","XY")]]

#get the centroid of the new cell
siteXspp_new<-cbind(siteXsppSUM,xyFromCell(a.proj,siteXsppSUM$cellnumber))

##I'm testing this out, taking this section out. 
#test <- colnames(siteXspp)
#test2 <- gsub('[.]','_',test)
#setnames(siteXspp,colnames(siteXspp),test2)


substructure<-siteXspp_new[,!colnames(siteXspp_new) %in% c("cellnumber","x","y"),with=F]
richness<-apply(substructure,1,sum)

a.proj[siteXspp_new$cellnumber]<-richness

plot(a.proj)

#write table to file
write.csv(siteXspp_new,paste(droppath,"Dimensions Data/siteXspp1degree.csv",sep="/"))
