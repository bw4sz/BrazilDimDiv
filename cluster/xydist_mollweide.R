####Find XY coordinates to all cells on the earth's surface

droppath<-"/work/02443/bw4sz/GlobalMammals/"

setwd(droppath)

#read in packages and source

require(reshape2)
require(maptools)
require(raster)
require(sp)

#read in xy data with original rows
xytab<-read.csv("Output/xytable.csv")[,-1]
xytab<-xytab[with(xytab, order(V1)), ] ## CP it has to be sorted by V1 to later find coordinates of each cell
rownames(xytab)<-xytab$V1

#drop the id row
f<-xytab[,colnames(xytab) %in% c("x","y")]

#transform into raster
f<-SpatialPoints(f)

crs(f)<-"+proj=moll +ellps=WGS84"
proj.f<-spTransform(f,crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "))

#calculate distances between points
xymat<-pointDistance(coordinates(proj.f),lonlat=TRUE)
diag(xymat)<-NA

#melt the distance matrix
xymelt<-melt(xymat)

head(xymelt)

xyout<-xymelt[!is.na(xymelt$value),]

#name the columns
colnames(xyout)<-c("To.OriginalRow","From.OriginalRow","km") #CP names changed this ID corresponds to v1 (original id)

print(dim(xyout))

print(head(xyout))

write.table(xyout,"xydist.txt")
