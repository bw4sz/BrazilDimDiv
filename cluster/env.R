suppressMessages(require(data.table))
suppressMessages(require(picante))

###Testing pbdMPI dist functions

require(pbdMPI,quietly=TRUE)

init()

#############################
######Source files from snoweye github

source_github <- function(u) {
  # load package
  suppressMessages(require(RCurl))
  
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

if(comm.rank()==0){
  env<-read.csv("Input/siteXsppXenv_1dgr.csv")
   
  env.table<-data.table(env[,colnames(env) %in% c("x","y",paste("bio",1:19,sep=""))])
  head(colnames(env))
  
  env.pca<-prcomp(env.table[,-c(1:2),with=F])$x

} else {

env.pca<-NULL
}

X.gbd <- comm.as.gbd(env.pca)

dist.X.common <- comm.dist(X.gbd)
dist.X.gbd <- comm.dist(X.gbd, return.type = "gbd")

### Verify 2.

comm.print(dist.X.gbd[1:10,])

print(dim(dist.X.gbd))

comm.write.table(dist.X.gbd,"Output/EnvData.txt")

finalize()
