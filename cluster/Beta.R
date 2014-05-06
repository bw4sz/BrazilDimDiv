#Taxonomic, Phylogenetic and Functional Diversity in Global Mammals

#Ben Gregory Weinstein, corresponding author - alll code below was writen by BGW

##################

require(pbdMPI,quietly=TRUE,warn.conflicts=FALSE)

init()


#Require libraries
suppressMessages(require(reshape,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(picante,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(ggplot2,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(parallel,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(foreach,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(GGally,quietly=TRUE,warn.conflicts=FALSE))
suppressMessages(require(data.table,quietly=TRUE,warn.conflicts=FALSE))

#Everyone say hello
comm.print(comm.rank(), all.rank = TRUE)

##If running locally set droppath

droppath<-"/home1/02443/bw4sz/GlobalMammals/"

setwd(droppath)

###Define Source Functions
suppressMessages(source("Input/BrazilSourceFunctions.R"))

comm.print(paste("Total nodes is:", comm.size()))

############ Read in Data

#set rank
.rank<-comm.rank()+1

#rank to filename
fil<-paste("Data/",.rank,"Rank.RData",sep="")

#read in portion of the data
load(fil)

timeF<-system.time(beta_out<-betaPar.scatter(comm.df,Index_Space,tcellbr_sub,traitm))

#Return timing argument to console
print(timeF)

#match to the xytable

# create all combinations of 2
# return as a data.table with these as columns `V1` and `V2`
# then count the numbers in each group
combs<-xytable[, rbindlist(combn(V1,2, 
                      FUN = function(x) as.data.table(as.list(x)), simplify = F))]
setkey(combs,V1)
setkey(xytable,V1)

#combine both ways
combsV1<-xytable[combs]
setkey(combsV1,V2)
combsV2<-xytable[combsV1]

#set names
setnames(combsV2,colnames(combsV2),c("To.OriginalRow","To.xcord","To.ycord","To","From.OriginalRow","From.xcord","From.ycord","From"))

#Bring in pairwise betadiveristy data
print(head(beta_out))

#set keys to match
setkey(beta_out,To,From)

#set to character
combsV2[,To:=as.character(To)]
combsV2[,From:=as.character(From)]

setkey(combsV2,To,From)

#create combo columns
combR<-apply(beta_out,1,function(x){paste(sort(c(as.numeric(x[1]),as.numeric(x[2]))),collapse="_")})


combRxy<-apply(combsV2,1,function(x){paste(sort(c(as.numeric(x[4]),as.numeric(x[8]))),collapse="_")})

beta_out[,combo:=combR,]
combsV2[,combo:=combRxy]

setkey(beta_out,combo)
setkey(combsV2,combo)

head(mergeT<-merge(combsV2,beta_out))

print(dim(mergeT))

#try writing from all
comm.write.table(mergeT,"Output/FinalData.csv",row.names=FALSE)

#remove data file
file.remove(fil)

finalize()

###############################
#Data Generation Complete
###############################