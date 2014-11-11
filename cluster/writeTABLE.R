

require(pbdMPI)

init()

suppressMessages(require(data.table))
suppressMessages(require(picante))
require(reshape2)
require(parallel)

###############
#setwd on cluser
droppath<-"/work/02443/bw4sz/GlobalMammals/"

setwd(droppath)

#I'm having trouble writing to file

load("Output/BetaDistEnv.RData")

ls()

print(gc())

tables()


a<-splitIndices(nrow(BetaEnvDist),comm.size())

towrite<-BetaEnvDist[a[[comm.rank()+1]],]

comm.write.table(towrite,"Output/BetaEnvDist.txt")

finalize()
