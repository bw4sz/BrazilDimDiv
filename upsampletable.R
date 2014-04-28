###Create a siteXspp table at a new resolution 

require(reshape)
require(foreach)
require(doSNOW)
require(raster)

##read in data
droppath<-"C:/Users/Jorge/Dropbox"

siteXspp<-read.csv(paste(droppath,"Dimensions Data/siteXspp.csv",sep="/"))
test <- colnames(siteXspp)
test2 <- gsub('[.]','_',test)
colnames(siteXspp) <- test2

#create a molleweide raster with the desired resolution
a<-raster()
res(a)<-1

#project raster
a.proj<-projectRaster(a,crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")

#Get XY coordinates
XY<-colsplit(siteXspp$X," ",c("X","Y"))

#get the cell number from the projected raster
cellnumber<-cellFromXY(a.proj,XY)

#create an index of XY to cell number
sp.index<-data.frame(XY=siteXspp$X,cellnumber)

#get a list of unique cell numbers
uC<-unique(cellnumber)

#Combine columns with same cells
#walk through each cell
cl<-makeCluster(5,"SOCK")
registerDoSNOW(cl)
newsiteXspp<-foreach(x=1:length(uC),.packages="raster") %dopar% {
  target<-uC[x]
  
  #get all the XY coords that fall in that cell
  targetXY<-sp.index[sp.index$cellnumber %in% target,"XY"]
  
  #get all the rows in the Xy targets
  tocombine<-siteXspp[siteXspp$X %in% targetXY,]
  
  #merge target rows
  combrow<-(apply(tocombine[,-1],2,sum) > 0)*1
  
  #get the centroid of the new cell
  newC<-xyFromCell(a.proj,target)
  out<-data.frame(cbind(XY=newC,t(as.matrix(combrow))))
  return(out)
}
stopCluster(cl)

#Bind together
newsiteXsppdf<-rbind.fill(newsiteXspp)

#visualize richness to check accuracy
cellnumberR<-cellFromXY(a.proj,cbind(newsiteXsppdf$x,y=newsiteXsppdf$y))

richness<-apply(newsiteXsppdf[,-c(1,2)],1,sum)

a.proj[cellnumberR]<-richness

plot(a.proj)

#write table to file
write.csv(newsiteXsppdf,paste(droppath,"Dimensions Data/siteXspp1degree.csv"))
