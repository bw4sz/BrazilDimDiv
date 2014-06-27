####test for geographic and environmental distances on combinations of high/low beta tax, phylo and traits

require(data.table)
require(multcomp)

droppath<-"C:/Users/Caterina/Documents/R/1312_beta/BrazilDimDiv"
setwd(droppath)

load(file="BetaDistEnv.RData")

hist(Beta$Sorenson) 
hist(Beta$envdist)
hist(Beta$MNTD) 
hist(Beta$BetaSim) 

setkey(Beta,Sorenson)

####select only cells without 0 and 1 beta tax values ??
Beta<-Beta[Sorenson!=0 & Sorenson!=1,]


###select 5% higher/lower beta tax
l5tax<-Beta[1:(nrow(Beta)*0.05),]
h5tax<-Beta[(nrow(Beta)-(nrow(Beta)*0.05)+1):nrow(Beta),]

summary(l5tax)
summary(h5tax)


####residuals of beta phylo~tax and trait~tax LOW TAX
##if data is too big try with: https://github.com/austian/bigExplore/blob/master/R/bigResiduals.R

#to check the shape of the relationship
library(MASS)
plot(gam(BetaSim~s(MNTD,k=4),data=l5tax))

#####Low tax residuals from phylo and traits
#phylo
resphy_ltax<-as.data.table(residuals(glm(BetaSim~Sorenson,data=l5tax)))
#plot(glm(BetaSim~Sorenson,data=l5tax))
setnames(resphy_ltax,"lrbphy")
hist(resphy_ltax$lrbphy)

#trait
restrait_ltax<-as.data.table(residuals(glm(MNTD~Sorenson,data=l5tax)))
#plot(glm(MNTD~Sorenson,data=l5tax))
setnames(restrait_ltax,"lrbtrait")
hist(restrait_ltax$lrbtrait)

allres_ltax<-cbind(l5tax,resphy_ltax,restrait_ltax)


####same for HIGH TAX     
#phylo
resphy_htax<-data.frame(residuals(glm(BetaSim~Sorenson,data=h5tax)))
#plot(glm(BetaSim~Sorenson,data=h5tax))
setnames(resphy_htax,"hrbphy")
hist(resphy_htax$hrbphy)

#trait
restrait_htax<-data.frame(residuals(glm(MNTD~Sorenson,data=h5tax)))
#plot(glm(BetaSim~MNTD,data=h5tax))
setnames(restrait_htax,"hrbtrait")
hist(restrait_htax$hrbtrait)

allres_htax<-cbind(h5tax,resphy_htax,restrait_htax)



###identify all combinations of phylo/trait: high/high high/low low/high low/low for LOW TAX
newcolname = "category"
allres_ltax[,newcolname:=NA,with=FALSE]
allres_ltax$category <-as.factor(ifelse(allres_ltax$lrbphy<0 & allres_ltax$lrbtrait<0,"lplt", 
                                        ifelse(allres_ltax$lrbphy>0 & allres_ltax$lrbtrait<0,"hplt",
                                               ifelse(allres_ltax$lrbphy<0 & allres_ltax$lrbtrait>0,"lpht",     
                                                      ifelse(allres_ltax$lrbphy>0 & allres_ltax$lrbtrait>0,"hpht","NA"))))) 


###identify all combinations of phylo/trait: high/high high/low low/high low/low for HIGH TAX
newcolname = "category"
allres_htax[,newcolname:=NA,with=FALSE]
allres_htax$category <-as.factor(ifelse(allres_htax$hrbphy<0 & allres_htax$hrbtrait<0,"lplt", 
                                        ifelse(allres_htax$hrbphy>0 & allres_htax$hrbtrait<0,"hplt",
                                               ifelse(allres_htax$hrbphy<0 & allres_htax$hrbtrait>0,"lpht",     
                                                      ifelse(allres_htax$hrbphy>0 & allres_htax$hrbtrait>0,"hpht","NA")))))

#####PLOT RESIDUALS
plot(allres_ltax$lrbphy,allres_ltax$lrbtrait) 
abline(v=0,h=0)
plot(allres_htax$hrbphy,allres_htax$hrbtrait)  
abline(v=0,h=0)

####ENVIRONMENT AND DISTANCE
##test for differences geographic distance between low tax and high tax
newcolname = "taxcat"
allres_ltax[,newcolname:="ltax",with=FALSE]
allres_htax[,newcolname:="htax",with=FALSE]
l<-list(allres_ltax,allres_htax)
allres<-rbindlist(l)
## allres<-system.time(rbind(allres_ltax,allres_htax,use.names=FALSE)) #same time with my subset but might be fatser with bigger data
boxplot(allres$km~allres$taxcat)
mod0<-glm(km~taxcat,data=allres)
summary(mod0)


#effect of geographic distance on category
#boxplot(km~category,data=allres_ltax)
mod1<-glm(km~category,data=allres_ltax)
summary(mod1)


#effect of environment distance on category
#boxplot(envdist~category,data=allres_ltax)
mod2<-glm(envdist~km+category,data=allres_ltax)
summary(mod2)
tuk<-glht(mod2,linfct=mcp(category="Tukey"))
summary(tuk)
tuk.cld<-cld(tuk)
plot(tuk.cld)





