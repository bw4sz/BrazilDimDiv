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
setwd(droppath)

### Examples.
comm.set.seed(10, diff = TRUE)

#create some fake data
X.gbd<-matrix(ncol=10,nrow=10,data=rbinom(100,1,.5))

dist.X.common <- comm.dist(X.gbd)
dist.X.gbd <- comm.dist(X.gbd, return.type = "gbd")

### Verify 1.
dist.X <- dist(do.call("rbind", allgather(X.gbd)))
comm.print(all(dist.X == dist.X.common))

### Verify 2.
dist.X.df <- do.call("rbind", allgather(dist.X.gbd))
comm.print(all(dist.X == dist.X.df[, 3]))
comm.print(dist.X)
comm.print(dist.X.df)

### Finish.
finalize()

#
finalize()
