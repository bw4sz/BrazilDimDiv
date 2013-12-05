#Taxonomic, Phylogenetic and Functional Diversity in Brazilian Mammals
#Written by Ben Weinstein - Stony Brook University

######################Below are all the packages needed for this work, if any of them flag an error, just say install.package( "name of package")
require(vegan)
require(reshape)
require(maptools)
require(raster)
require(rgdal)
require(ape)
require(picante)
require(ggplot2)
require(gdistance)
library(doSNOW)
library(foreach)
require(stringr)
require(scales)

##################
#Set Paths
##################

#Set dropbox path
droppath<-"C:/Users/Ben/Dropbox/Dimensions/Data/Brazil"

#Set git path
gitpath<-"C:/Users/Ben/Documents/BrazilDimDiv/"

#sink output for overnight runs so we can see it tomorrow
#sink("")
#load data if desired

#load(paste(droppath,"EnvCompare.RData",sep=""))

###Define Source Functions
source(paste(gitpath,"BrazilSourceFunctions.R",sep=""))

###########################################
#Read in data
###########################################

#Read in species matrix
siteXspp <- read.csv(paste(drop"siteXspp_Oct20_2011.csv", row.names=1)

#Get entire species list
splist<-colnames(siteXspp)

#Read in phylogeny
tree<-read.nexus(paste(droppath,"New_Mammal_Timetree_IUCN_2Feb2013.nwk",sep=""))

#bring in traits
mon <- read.csv("C:\\Users\\Jorge\\Dropbox\\Lab paper 1 Predicted vs observed assemblages\\MorphologyShort.csv",na.strings="9999")

#make sure to match the tip labels and the siteXspp matrix
table(tree$tiplabels %in% colnames(siteXspp))
#######################################################
#Compute Environmental Dissimilarity
#######################################################
#Step 2 Bring in Env Info and Dissimilarity 
##############################################################
comm<-siteXspp

#Extract env information from each locality

#Either import a file showing the XY points of your localities
#loc<-readShapePoints("Locs_proj.shp")
#head(loc)
mean(edist<-pointDistance(loc,longlat=TRUE)/1000),na.rm=TRUE)
mean(edist)

#Or if the rownames of the siteXspp matrix are the cells that match the env layer
loc<-xyFromCell(rownames(siteXspp))

#Put all raster layers into one folder, point towards that folder.
clim.var <- stack(paste(droppath,'Env_Layers/BioClimBrazil.grd',sep="/")) ##open bioclim variables

#Aggregate if desired to the cell size of the siteXspp matrix, they need to line up!
#clim.var.1.5dgr <- aggregate(clim.var, fact=3, fun=mean) ### change it to 1.5 degree cells
plot(clim.var.1.5dgr[[9]])

remote <- stack(list.files(paste(droppath,"Env_Layers",sep="/"),pattern='.tif',full.names=TRUE)) ## create a stack of remote sensing variables
projection(remote) <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0' ##define projection

#aggregate to level of the climate layers
remote.1.5dgr <- aggregate(remote, fact=6, fun=mean) ### scale up to 1.5 degree cells
#get the identical extent of the climate layers
e <- extent(clim.var)
remote.sens <- crop(remote.1.5dgr, e)
remote.sens <- resample(remote.sens, clim.var.1.5dgr) ###make sure remote sense aligns with bioclim

plot(remote.sens[[1]])
env.list <- stack(clim.var.1.5dgr, remote.sens) ###join all environmental variables in one stack

#Name the environmental layers
names(env.list)

#Extract raster info to points
extract.list<-lapply(env.list, function(x) extract(x,loc))

Envtable<-data.frame(loc@data,as.data.frame(extract.list))

#remember to divide temp by 10
Envtable$AnnualTemp<-Envtable$AnnualTemp/10

#check for any NA's
IDerror<-Envtable[!complete.cases(Envtable),]

#Remove those NA's if needed
#Envtable<-Envtable[!complete.cases(Envtable),]

#WRite to file
write.csv(Envtable,paste(gitpath,"Envtable.csv",sep=","))

#Create distance matrices of each Environmental variable, store it in a list. 
Evar<-lapply(Envtable[7:12],function(x)(dist(x, method='maximum', diag=TRUE, upper=TRUE,)))

#Create euclidean distance matrix from the long/lat coordinates
longlat<-as.matrix(cbind(Envtable$LongDecDeg,Envtable$LatDecDeg))
rownames(longlat)<-Envtable$ID_Comm
edist<-rdist.earth(longlat,longlat,miles=FALSE)
Evar<-lappend(Evar,edist)

#Bring in Cost Path from a R gdistance script
#############iF you are not changing ANY of the above steps uncomment the next 2 lines, and SKIP this section of data
#CostPathMatrix<-read.csv("CostPathCost.csv", row.names=1)
#rownames(CostPathMatrix)<-Envtable$ID_Comm
#colnames(CostPathMatrix)<-Envtable$ID_Comm

#Import Friction layer
elevr<-raster(paste(droppath,"",sep=""))

#Try different elevation layer
elev.test<-crop(elevr,extent(loc)*1.2)

#Desired res? this needs to be checked carefully, as res decreases huge computation increases
elev.test<-aggregate(elev.test,100)
res(elev.test)

#Find shortest cost path
cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)
costPath.list<-foreach(x = 1:length(loc),.packages=c("raster","gdistance")) %dopar% {
  
  #pick the original site. 
  orig<-loc[x,]
  #What elevation is the origin
  elev_origin<-extract(elev.test,orig)[[1]]
  
  #Get the difference between the origin elevation and every cell in the raster
  elev_diff<-abs(elev_origin-elev.test)
  
  #create a the transition matrix where permeability is high where elev difference is low
  trans_diff<-transition(elev_diff,function(x) 1/mean(x),8)
  
  #Diagonal Cell Correction, diagonals are father away than adjacent cells. 
  slope <- geoCorrection(trans_diff)
  
  #Remember this cost surface is only valid for this site as the origen, ie. we want to create just this column in the cost distance matrix
  #Cost Distance
  cdist<-costDistance(slope,orig,loc)
  #labelling is key here.
  
  return(list(cdist))}
stopCluster(cl)

#reshape cost path for each locality pair
m.costlist<-melt(costPath.list)
CostPathMatrix<-cast(m.costlist,X2~L1)[,-1]
rownames(CostPathMatrix)<-loc$ID_Comm
colnames(CostPathMatrix)<-loc$ID_Comm

#Write to file
write.csv(CostPathMatrix,paste(gitpath,"CostPathCost.csv",sep="")

Evar<-lappend(Evar,as.matrix(CostPathMatrix))
names(Evar)<-c(names(Envtable[7:ncol(Envtable)]),"Euclid", "CostPathCost")

#Write Environmental Layers Matrices to File
setwd(paste(gitpath,"Results",sep="/"))
          
for(x in 1:length(Evar))
  write.csv(as.matrix(Evar[[x]]),file=paste(names(Evar[x]),".csv",""))

#create giant data frame using the comparisons between all sites
Evar.m<-lapply(Evar,function(x) {
  y<-as.matrix(x)
  colnames(y)<-Envtable$ID_Comm
  rownames(y)<-Envtable$ID_Comm
  return(y)})

test.Evar<-melt(Evar.m)
compare.env<-cast(test.Evar,X1+X2~L1)
colnames(compare.env)<-c("To","From",names(compare.env[-c(1,2)]))

#save final matrix
write.csv(compare.env,"EnvCompare.csv")

#save workspace
save.image(paste(droppath,"EnvCompare.RData",sep=""))
