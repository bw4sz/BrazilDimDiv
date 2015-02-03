require(data.table)

#setwd
droppath<-"/work/02443/bw4sz/GlobalMammals/"

print("Beta")
print(dat)

###Step 1## Match unique ID to richness table
# ############################################################################
siteXrich<-fread("Input/siteXrich.csv")

print("siteXrich")
print(siteXrich)

# create all combinations of 2
# return as a data.table with these as columns `V1` and `V2`
# then count the numbers in each group

system.time(combs<-as.data.table(t(combn(siteXrich$V1,2, simplify = T))))

print(combs)

print(gc())

#Merge richness info, set keys
setkey(combs,V1)
setkey(siteXrich,V1)

#combine both ways
combsV1<-siteXrich[combs]
setkey(combsV1,V2)
combsV2<-siteXrich[combsV1]

#set names
setnames(combsV2,colnames(combsV2),c("From.OriginalRow","From.Richness","To.OriginalRow","To.Richness"))

print("richtable")
print(combsV2)
rm(combsV1)
rm(combs)
rm(siteXrich)
gc()


###Step 2## Match with Randomization table
#########################################################################
ran<-fread("BrazilDimDiv4deg/Randomization4dg.txt")
setnames(ran,names(ran),c("Phylo","Tax","Trait","Quantile","From.Richness","To.Richness"))
print("randomization table")
print(ran)

#Reshape randomization table to have a unique to.rich/from.rich id
ran2<-reshape(ran, idvar=c("From.Richness","To.Richness"), timevar='Quantile', direction='wide')

#set key
setkey(combsV2,From.Richness,To.Richness)
setkey(ran2,From.Richness,To.Richness)

#first merge
mergeL<-merge(combsV2,ran2)

#Switch to and from column name 
setnames(ran2,colnames(ran2),c("To.Richness","From.Richness","Phylo.Lower","Tax.Lower","Trait.Lower","Phylo.Upper","Tax.Upper","Trait.Upper"))
setkey(ran2,From.Richness,To.Richness)

#merge again
mergeR<-merge(combsV2,ran2)

rm(combsV2)
rm(ran)
rm(ran2)
print(gc())

datrich<-rbindlist(list(mergeR,mergeL))
setkey(datrich,From.OriginalRow,To.OriginalRow)
datrich<-unique(datrich)

print("datrich")
print(datrich)

rm(mergeR)
rm(mergeL)
gc()


###Step 3## Merge with BetaEnvDist table
#########################################################################

# #Betadiversity data
load(BetaEnvDist.RData)
print(head(BetaEnvDist))

setkey(BetaEnvDist,From.OriginalRow,To.OriginalRow)

a<-merge(BetaEnvDist,datrich)

print(a)

#reverse names
setnames(datrich,colnames(datrich),c("From.Richness","To.Richness","To.OriginalRow","From.OriginalRow",
                                     "Phylo.Lower","Tax.Lower","Trait.Lower",     
                                     "Phylo.Upper","Tax.Upper","Trait.Upper"))

#set keys
setkey(datrich,From.OriginalRow,To.OriginalRow)

b<-merge(BetaEnvDist,datrich)

d<-rbindlist(list(a,b))

print("Final Table")
print(d)

rm(a)
rm(b)
rm(datrich)

print(gc())

BetaEnvDist<-unique(d)

print(BetaEnvDist)

rm(d)
gc()

save.image("Output/BetaEnvDistN.RData")

write.table(BetaEnvDist,"Output/BetaDistEnvN.txt",row.names=FALSE)

#names(BetaEnvDist)
#[1] "From.OriginalRow" "To.OriginalRow"   "envdist"          "km"               "From.xcord"       "From.ycord"       "To.xcord"        
#[8] "To.ycord"         "Phylo"            "Tax"              "Trait"            "From.Richness"    "To.Richness"      "Phylo.Lower"     
#[15] "Tax.Lower"        "Trait.Lower"      "Phylo.Upper"      "Tax.Upper"        "Trait.Upper"     
