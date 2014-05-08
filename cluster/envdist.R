require(data.table)
require(picante)

###Testing pbdMPI dist functions

require(pbdMPI,quietly=TRUE)

init()

#############################
######Source files from snoweye github

source_github <- function(u) {
  # load package
  require(RCurl)
  
  # read script lines from website and evaluate
  script <- getURL(u, ssl.verifypeer = FALSE)
  eval(parse(text = script),envir=.GlobalEnv)
}  

#source allpairs.R
source_github("https://raw.githubusercontent.com/snoweye/pbdMPI/master/R/comm_allpairs.r")

#source comm.dist R
source_github("https://raw.githubusercontent.com/snoweye/pbdMPI/master/R/comm_dist.r")

#source comm as gbd
source_github("https://raw.githubusercontent.com/snoweye/pbdMPI/master/R/comm_as_gbd.r")

#looks like we also need the comm.allcommon tools from this file
source_github("https://raw.githubusercontent.com/snoweye/pbdMPI/master/R/comm_tool.r")

###############
#setwd on cluser
droppath<-"/home1/02443/bw4sz/GlobalMammals/"

##If running locally set dropbox path
#droppath<-"C:/Users/Jorge/Documents/BrazilDimDiv/cluster"
setwd(droppath)

### Initial.
library(pbdMPI, quietly = TRUE)
init()



dist(env.pca$x[1:10,])

if(comm.rank()==0){
  env<-read.csv("Input/siteXsppXenv_1dgr.csv",nrows=1000)
  
  tree<-read.tree("/Input//Sep19_InterpolatedMammals_ResolvedPolytomies.nwk")
  
  env.table<-data.table(env[,colnames(env) %in% c("x","y",paste("bio",1:19,sep=""))])
  head(colnames(env))
  
  env.pca<-prcomp(env.table[,-c(1:2),with=F])
}

X.gbd <- comm.as.gbd(env.pca$x)

dist.X.common <- comm.dist(X.gbd)
dist.X.gbd <- comm.dist(X.gbd, return.type = "gbd")

### Verify 2.
dist.X.df <- do.call("rbind", allgather(dist.X.gbd))

comm.print(dist.X.df[1:10,])

comm.write.table(mergeT,"Output/EnvData.csv")
