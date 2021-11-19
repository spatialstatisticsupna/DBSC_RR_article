rm(list=ls())
library(bigDM)
library(INLA)
library(sf)
library(RColorBrewer)
library(tmap)


#########################################
## Load the data and cartography files ##
#########################################
load("../data/England_MSOA.Rdata")
head(carto)

S <- nrow(carto)

Qs <- Diagonal(S,colSums(W))-W
Qs.Leroux <- Diagonal(S)-Qs


#########################
## Auxiliary functions ##
#########################
if(!file.exists("figures")) {
  dir.create(file.path(getwd(), "figures"))
}

DSS.score <- function(model){
  
  S <- nrow(model$.args$data)
  
  O.hat <- numeric(S)
  O.sd <- numeric(S)
  
  for(i in 1:S){
    E.r <- inla.emarginal(function(x) x, model$marginals.fitted.values[[i]])
    Var.r <- inla.emarginal(function(x) x^2, model$marginals.fitted.values[[i]])-E.r^2
    O.hat[i] <- model$.args$data$E[i]*E.r
    O.sd[i] <- sqrt(model$.args$data$E[i]*E.r+model$.args$data$E[i]^2*Var.r)
  }
  
  DSS <- ((model$.args$data$O-O.hat)/O.sd)^2+2*log(O.sd)
  return(sum(DSS))
}

DIC <- function(x){
  data.frame(mean.deviance=x$dic$mean.deviance, ## posterior mean deviance
             pD=x$dic$p.eff,                    ## effective number of parameters
             DIC=x$dic$dic,                     ## Deviance Information Criterion
             WAIC=x$waic$waic,                  ## Watanabe-Akaike information criterion
             LS=-sum(log(x$cpo$cpo)),           ## Logarithmic Score (see inla.cpo function)
             DSS=DSS.score(x))                  ## Dawid-Sebastiani Score
}


###########################################################################################################
## Load all the models fitted with INLA, which can be downloaded from:                                   ##
##  - https://emi-sstcdapp.unavarra.es/England_MSOA/INLA_models/England_MSOA_SpatialModels_iCAR.Rdata    ##
##  - https://emi-sstcdapp.unavarra.es/England_MSOA/INLA_models/England_MSOA_SpatialModels_LCAR.Rdata    ##
##  - https://emi-sstcdapp.unavarra.es/England_MSOA/INLA_models/England_MSOA_SpatialModels_BYM.Rdata     ##
##  - https://emi-sstcdapp.unavarra.es/England_MSOA/INLA_models/England_MSOA_SpatialModels_BYM2.Rdata    ##
##  - https://emi-sstcdapp.unavarra.es/England_MSOA/INLA_models/England_MSOA_ClusterModels_iCAR.Rdata    ##
##  - https://emi-sstcdapp.unavarra.es/England_MSOA/INLA_models/England_MSOA_ClusterModels_LCAR.Rdata    ##
##  - https://emi-sstcdapp.unavarra.es/England_MSOA/INLA_models/England_MSOA_ClusterModels_BYM.Rdata     ##
##  - https://emi-sstcdapp.unavarra.es/England_MSOA/INLA_models/England_MSOA_ClusterModels_BYM2.Rdata    ##
##  - https://emi-sstcdapp.unavarra.es/England_MSOA/INLA_models/England_MSOA_ClusterModels_BYM2_RR.Rdata ##
###########################################################################################################
load("England_MSOA_SpatialModels_iCAR.Rdata")
load("England_MSOA_SpatialModels_LCAR.Rdata")
load("England_MSOA_SpatialModels_BYM.Rdata")
load("England_MSOA_SpatialModels_BYM2.Rdata")
load("England_MSOA_ClusterModels_BYM2.Rdata")
load("England_MSOA_ClusterModels_iCAR.Rdata")
load("England_MSOA_ClusterModels_LCAR.Rdata")
load("England_MSOA_ClusterModels_BYM.Rdata")
load("England_MSOA_ClusterModels_BYM2.Rdata")
load("England_MSOA_ClusterModels_BYM2_RR.Rdata")


###########################################################################################
## Fig 1: Posterior marginal distributions of the regression coefficients (BYM2+C model) ##
###########################################################################################
MODELS <- list(GLM,BYM2.C$`l=1`,BYM2.C.RR)
names(MODELS) <- c("GLM","BYM2+C","BYM2+C+RR")

title <- c("Lieberson index (ISOL)","Nursing homes (NH)","Health Deprivation and Disability (HDD)","Air quality (AIRQ)")
pos <- c("topleft","topright","topleft","topleft")

graphics.off()
pdf("figures/Figure1.pdf", width=8, height=6)
par(mfrow=c(2,2), pty="m")

for(i in 1:5){
  plot(inla.smarginal(MODELS$GLM$marginals.fixed[[i+1]]), type="l", main=title[i], xlab="", ylab="", xlim=c(0,0.25))
  lines(inla.smarginal(MODELS$`BYM2+C`$marginals.fixed[[i+1]]), col="red")
  lines(inla.smarginal(MODELS$`BYM2+C+RR`$marginals.fixed[[i+1]]), col="blue")
  legend(pos[i], legend=names(MODELS), col=c("black","red","blue"), lwd=1, bty="n")
}
dev.off()


###########################################################################
## Fig 2: Maps with the posterior median estimates of relative risks and ##
##        posterior exceedence probabilities obtained with BYM2+C model  ##
###########################################################################
carto$risk.BYM2.C <- BYM2.C.RR$summary.fitted.values$`0.5quant`
carto$prob.BYM2.C <- 1-BYM2.C.RR$summary.fitted.values$`1 cdf`

## Map of posterior median estimates of relative risks ##
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.67,0.77,0.91,1,1.10,1.30,1.50,Inf)

BYM2.C.risk.map <- tm_shape(carto) +
  tm_polygons(col="risk.BYM2.C", id="CODE", palette=paleta, border.alpha=0, title="RR",
              leyend.show=T, legend.reverse=T, style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Posterior median estimates of relative risks", main.title.position="left",
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01))

tmap_save(BYM2.C.risk.map, file="figures/Map_EstimatedRisks_BYM2+C+RR.pdf")


## Map of posterior exceedence probabilities ##
paleta <- brewer.pal(7,"Blues")[-c(1,5)]
values <- c(0,0.1,0.2,0.8,0.9,1)

BYM2.C.prob.map <- tm_shape(carto) +
  tm_polygons(col="prob.BYM2.C", id="CODE", palette=paleta, border.alpha=0, title="Probs",
              leyend.show=T, legend.reverse=T, style="fixed", breaks=values, interval.closure="left",
              labels=c("[0-0.1)","[0.1-0.2)","[0.2-0.8)","[0.8-0.9)","[0.9-1]")) +
  tm_layout(main.title="Posterior exceedence probabilities", main.title.position=0.1,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01))

tmap_save(BYM2.C.prob.map, file="figures/Map_ExceedenceProbs_BYM2+C+RR.pdf")


####################################################################################
## Online figure: Maps with the posterior median estimates of relative risks and  ##
##                posterior exceedence probabilities obtained with BYM2 models    ##
####################################################################################
carto$risk.BYM2 <- BYM2.RR$summary.fitted.values$`0.5quant`
carto$prob.BYM2 <- 1-BYM2.RR$summary.fitted.values$`1 cdf`

## Map of posterior median estimates of relative risks ##
paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.67,0.77,0.91,1,1.10,1.30,1.50,Inf)

BYM2.risk.map <- tm_shape(carto) +
  tm_polygons(col="risk.BYM2", id="CODE", palette=paleta, border.alpha=0, title="RR",
              leyend.show=T, legend.reverse=T, style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Posterior median estimates of relative risks", main.title.position="left",
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01))

tmap_save(BYM2.risk.map, file="figures/Map_EstimatedRisks_BYM2+RR.pdf")


## Map of posterior exceedence probabilities ##
paleta <- brewer.pal(7,"Blues")[-c(1,5)]
values <- c(0,0.1,0.2,0.8,0.9,1)

BYM2.prob.map <- tm_shape(carto) +
  tm_polygons(col="prob.BYM2", id="CODE", palette=paleta, border.alpha=0, title="Probs",
              leyend.show=T, legend.reverse=T, style="fixed", breaks=values, interval.closure="left",
              labels=c("[0-0.1)","[0.1-0.2)","[0.2-0.8)","[0.8-0.9)","[0.9-1]")) +
  tm_layout(main.title="Posterior exceedence probabilities", main.title.position=0.1,
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01))

tmap_save(BYM2.prob.map, file="figures/Map_ExceedenceProbs_BYM2+RR.pdf")


###########################################################################################
## Fig A.1: Posterior marginal distributions of the regression coefficients (iCAR model) ##
###########################################################################################
MODELS <- list(GLM,iCAR,iCAR.RR)
names(MODELS) <- c("GLM","iCAR","iCAR+RR")

title <- c("Lieberson index (ISOL)","Nursing homes (NH)","Health Deprivation and Disability (HDD)","Air quality (AIRQ)")
pos <- c("topleft","topright","topleft","topleft")

graphics.off()
pdf("figures/FigureA1.pdf", width=8, height=6)
par(mfrow=c(2,2), pty="m")

for(i in 1:5){
  plot(inla.smarginal(MODELS$GLM$marginals.fixed[[i+1]]), type="l", main=title[i], xlab="", ylab="", xlim=c(0,0.25))
  lines(inla.smarginal(MODELS$`iCAR`$marginals.fixed[[i+1]]), col="red")
  lines(inla.smarginal(MODELS$`iCAR+RR`$marginals.fixed[[i+1]]), col="blue")
  legend(pos[i], legend=names(MODELS), col=c("black","red","blue"), lwd=1, bty="n")
}
dev.off()


###########################################################################################
## Fig A.2: Posterior marginal distributions of the regression coefficients (LCAR model) ##
###########################################################################################
MODELS <- list(GLM,LCAR,LCAR.RR)
names(MODELS) <- c("GLM","LCAR","LCAR+RR")

title <- c("Lieberson index (ISOL)","Nursing homes (NH)","Health Deprivation and Disability (HDD)","Air quality (AIRQ)")
pos <- c("topleft","topright","topleft","topleft")

graphics.off()
pdf("figures/FigureA2.pdf", width=8, height=6)
par(mfrow=c(2,2), pty="m")

for(i in 1:5){
  plot(inla.smarginal(MODELS$GLM$marginals.fixed[[i+1]]), type="l", main=title[i], xlab="", ylab="", xlim=c(0,0.25))
  lines(inla.smarginal(MODELS$`LCAR`$marginals.fixed[[i+1]]), col="red")
  lines(inla.smarginal(MODELS$`LCAR+RR`$marginals.fixed[[i+1]]), col="blue")
  legend(pos[i], legend=names(MODELS), col=c("black","red","blue"), lwd=1, bty="n")
}
dev.off()


##########################################################################################
## Fig A.3: Posterior marginal distributions of the regression coefficients (BYM model) ##
##########################################################################################
MODELS <- list(GLM,BYM,BYM.RR)
names(MODELS) <- c("GLM","BYM","BYM+RR")

title <- c("Lieberson index (ISOL)","Nursing homes (NH)","Health Deprivation and Disability (HDD)","Air quality (AIRQ)")
pos <- c("topleft","topright","topleft","topleft")

graphics.off()
pdf("figures/FigureA3.pdf", width=8, height=6)
par(mfrow=c(2,2), pty="m")

for(i in 1:5){
  plot(inla.smarginal(MODELS$GLM$marginals.fixed[[i+1]]), type="l", main=title[i], xlab="", ylab="", xlim=c(0,0.25))
  lines(inla.smarginal(MODELS$`BYM`$marginals.fixed[[i+1]]), col="red")
  lines(inla.smarginal(MODELS$`BYM+RR`$marginals.fixed[[i+1]]), col="blue")
  legend(pos[i], legend=names(MODELS), col=c("black","red","blue"), lwd=1, bty="n")
}
dev.off()


###########################################################################################
## Fig A.4: Posterior marginal distributions of the regression coefficients (BYM2 model) ##
###########################################################################################
MODELS <- list(GLM,BYM2,BYM2.RR)
names(MODELS) <- c("GLM","BYM2","BYM2+RR")

title <- c("Lieberson index (ISOL)","Nursing homes (NH)","Health Deprivation and Disability (HDD)","Air quality (AIRQ)")
pos <- c("topleft","topright","topleft","topleft")

graphics.off()
pdf("figures/FigureA4.pdf", width=8, height=6)
par(mfrow=c(2,2), pty="m")

for(i in 1:5){
  plot(inla.smarginal(MODELS$GLM$marginals.fixed[[i+1]]), type="l", main=title[i], xlab="", ylab="", xlim=c(0,0.25))
  lines(inla.smarginal(MODELS$`BYM2`$marginals.fixed[[i+1]]), col="red")
  lines(inla.smarginal(MODELS$`BYM2+RR`$marginals.fixed[[i+1]]), col="blue")
  legend(pos[i], legend=names(MODELS), col=c("black","red","blue"), lwd=1, bty="n")
}
dev.off()


###########################################################################
## Table 1: Descriptive statistics of predictor variables (risk factors) ##
###########################################################################
Table1 <- data.frame(ISOL=c(mean(carto$ISOL),sd(carto$ISOL),min(carto$ISOL),max(carto$ISOL)),
                     NH=c(mean(carto$NH),sd(carto$NH),min(carto$NH),max(carto$NH)),
                     HDD=c(mean(carto$HDD),sd(carto$HDD),min(carto$HDD),max(carto$HDD)),
                     AIRQ=c(mean(carto$AIRQ),sd(carto$AIRQ),min(carto$AIRQ),max(carto$AIRQ)))
rownames(Table1) <- c("Average","Standard deviation","Minimum","Maximum")

round(Table1,3)


##################################################################
## Table 2: Model selection criteria fo models fitted with INLA ##
##################################################################
MODELS <- list(GLM,iCAR,iCAR.RR,LCAR,LCAR.RR,BYM,BYM.RR,BYM2,BYM2.RR)
names(MODELS) <- c("GLM","iCAR","iCAR+RR","LCAR","LCAR+RR","BYM","BYM+RR","BYM2","BYM2+RR")

Table2 <- do.call(rbind,lapply(MODELS, DIC))
round(Table2,1)


########################################################################
## Table 3: Posterior mean, posterior standard deviation, and         ##
##          95% credible intervals of the regression coefficients     ##
##          for GLM, ICAR, LCAR, BYM and BYM2 models fitted with INLA ##
########################################################################
MODELS <- list(GLM,iCAR,LCAR,BYM,BYM2,iCAR.RR,LCAR.RR,BYM.RR,BYM2.RR)
names(MODELS) <- c("GLM","iCAR","LCAR","BYM","BYM2","iCAR+RR","LCAR+RR","BYM+RR","BYM2+RR")

Table3 <- vector("list",5)
names(Table3) <- c("(Intercept)","ISOL","NH","HDD","AIRQ")

for(i in seq(length(Table3))){
  Table3[[i]] <- do.call(rbind, lapply(MODELS, function(x) round(x$summary.fixed[i,1:5],3)))
  rownames(Table3[[i]]) <- names(MODELS)
}

print(Table3[-1])


########################################################################
## Table 4: Posterior mean, posterior standard deviation, and         ##
##          95% credible intervals of the model hyperparameters       ##
##          for GLM, ICAR, LCAR, BYM and BYM2 models fitted with INLA ##
########################################################################
hyper.iCAR <- iCAR$summary.hyperpar[,1:5]
hyper.LCAR <- LCAR$summary.hyperpar[,1:5]
hyper.BYM <- BYM$summary.hyperpar[,1:5]
hyper.BYM2 <- BYM2$summary.hyperpar[,1:5]

Table4 <- rbind(hyper.iCAR,hyper.LCAR,hyper.BYM,hyper.BYM2)
rownames(Table4) <- c("iCAR - tau.s","LCAR - tau.s","LCAR - lambda.s",
                      "BYM - tau.u","BYM - tau.v","BYM2 - tau.s","BYM2 - lambda.s")

round(Table4,3)


##########################################################################
## Table 5: Model selection criteria for BYM2+C models fitted with INLA ##
##          considering different neighborhood orders (parameter l)     ##
##########################################################################
Table5 <- do.call(rbind,lapply(BYM2.C, DIC))
round(Table5,1)


######################################################################
## Table 6: Posterior mean, posterior standard deviation, and 95%   ##
##          credible intervals of the regression coefficients for   ##
##          GLM, BYM2+C and BYM2+C+RR models (l=1) fitted with INLA ##
######################################################################
MODELS <- list(GLM,BYM2.C$`l=1`,BYM2.C.RR)
names(MODELS) <- c("GLM","BYM2+C","BYM2+C+RR")

Table6 <- vector("list",5)
names(Table6) <- c("(Intercept)","ISOL","NH","HDD","AIRQ")

for(i in seq(length(Table6))){
  Table6[[i]] <- do.call(rbind, lapply(MODELS, function(x) round(x$summary.fixed[i,1:5],3)))
  rownames(Table6[[i]]) <- names(MODELS)
}

print(Table6[-1])


######################################################################################
## Table 7: Total MSOAs classified as extreme relative risk by Urban-Rural category ##
######################################################################################

## Compute average relative risks in neighbouring MSOAs ##
nb.risk <- function(i){
  
  pos <- as.numeric(which(W[i,]==1))
  w <- rep(1/length(pos),length(pos))
  
  post.mean <- model.summary.fitted.values[pos,"mean"]
  post.sd <- model.summary.fitted.values[pos,"sd"]
  post.cdf <- model.summary.fitted.values[pos,"1 cdf"]
  marginals <- model.marginals.fitted.values[pos]
  
  xx <- sort(unlist(lapply(marginals, function(x) x[,"x"])))
  at <- round(seq(1,length(xx),length.out=75))
  nb.marginals <- matrix(0, nrow=0, ncol=2, dimnames=list(NULL, c("x","y")))
  for(j in xx[at]){
    aux <- unlist(lapply(marginals, function(x) inla.dmarginal(j,x)))
    nb.marginals <- rbind(nb.marginals,c(j,sum(aux*w)))
  }
  
  nb.summary <- data.frame(sum(post.mean*w),
                           sqrt(sum((post.sd^2+post.mean^2-sum(post.mean*w)^2)*w)),
                           inla.qmarginal(0.025,nb.marginals),
                           inla.qmarginal(0.5,nb.marginals),
                           inla.qmarginal(0.975,nb.marginals),
                           inla.mmarginal(nb.marginals),
                           sum(post.cdf*w))
  colnames(nb.summary) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode","1 cdf")
  
  return(list(nb.summary=nb.summary,nb.marginals=nb.marginals))
}

## BYM2 model (this might take a while...) ##
ID <- as.list(1:S)
model.summary.fitted.values <- BYM2.RR$summary.fitted.values
model.marginals.fitted.values <- BYM2.RR$marginals.fitted.values

cl <- makeCluster(detectCores())
clusterExport(cl, varlist=c("model.summary.fitted.values","model.marginals.fitted.values","ID","W"), envir=environment())
clusterEvalQ(cl, library(INLA))
system.time(aux <- parLapply(cl,ID,nb.risk))
stopCluster(cl)

BYM2.nb <- vector("list",2)
names(BYM2.nb) <- c("nb.summary", "nb.marginals")

BYM2.nb$nb.summary <- do.call(rbind,lapply(aux, function(x) x$nb.summary))
BYM2.nb$nb.marginals <- lapply(aux, function(x) x$nb.marginals)
carto$BYM2.nb.prob1 <- 1-BYM2.nb$nb.summary$`1 cdf`


## BYM2+C model (this might take a while...) ##
ID <- as.list(1:S)
model.summary.fitted.values <- BYM2.C.RR$summary.fitted.values
model.marginals.fitted.values <- BYM2.C.RR$marginals.fitted.values

cl <- makeCluster(detectCores())
clusterExport(cl, varlist=c("model.summary.fitted.values","model.marginals.fitted.values","ID","W"), envir=environment())
clusterEvalQ(cl, library(INLA))
system.time(aux <- parLapply(cl,ID,nb.risk))
stopCluster(cl)

BYM2.C.nb <- vector("list",2)
names(BYM2.C.nb) <- c("nb.summary", "nb.marginals")

BYM2.C.nb$nb.summary <- do.call(rbind,lapply(aux, function(x) x$nb.summary))
BYM2.C.nb$nb.marginals <- lapply(aux, function(x) x$nb.marginals)
carto$BYM2.C.nb.prob1 <- 1-BYM2.C.nb$nb.summary$`1 cdf`
  

## Compute posterior exceedence probabilities of being greater than 1 or 1.5 ##
carto$BYM2.prob1 <- 1-BYM2.RR$summary.fitted.values$`1 cdf`
carto$BYM2.C.prob1 <- 1-BYM2.C.RR$summary.fitted.values$`1 cdf`

carto$BYM2.prob15 <- unlist(lapply(BYM2.RR$marginals.fitted.values, function(x) 1-inla.pmarginal(1.5, x)))
carto$BYM2.C.prob15 <- unlist(lapply(BYM2.C.RR$marginals.fitted.values, function(x) 1-inla.pmarginal(1.5, x)))


## Table 7 ##
Obs.SMR <- round(tapply(carto$O, carto$`Urban-Rural`, sum)/tapply(carto$E, carto$`Urban-Rural`, sum),2)
High.BYM2 <- tapply(carto$BYM2.prob15, carto$`Urban-Rural`, function(x) sum(x>0.9))
High.BYM2.C <- tapply(carto$BYM2.C.prob15, carto$`Urban-Rural`, function(x) sum(x>0.9))
Deaths.BYM2 <- tapply(carto[carto$BYM2.prob15>0.9,]$O, carto[carto$BYM2.prob15>0.9,]$`Urban-Rural`, sum)
Deaths.BYM2.C <- tapply(carto[carto$BYM2.C.prob15>0.9,]$O, carto[carto$BYM2.C.prob15>0.9,]$`Urban-Rural`, sum)
Overlapping.BYM2 <- table(carto[carto$BYM2.prob1>0.9 & carto$BYM2.nb.prob1>0.9,]$`Urban-Rural`)
Overlapping.BYM2.C <- table(carto[carto$BYM2.C.prob1>0.9 & carto$BYM2.C.nb.prob1>0.9,]$`Urban-Rural`)

Table7 <- cbind(Obs.SMR,High.BYM2,High.BYM2.C,Deaths.BYM2,Deaths.BYM2.C,Overlapping.BYM2,Overlapping.BYM2.C)
Table7[is.na(Table7)] <- 0
Table7 <- rbind(Table7,c(1,apply(Table7,2,sum)[-1]))
Table7


################################################################################
## Table 8: Total MSOAs classified as extreme relative risk by English region ##
################################################################################
Obs.SMR <- round(tapply(carto$O, carto$`Region`, sum)/tapply(carto$E, carto$`Region`, sum),2)
High.BYM2 <- tapply(carto$BYM2.prob15, carto$`Region`, function(x) sum(x>0.9))
High.BYM2.C <- tapply(carto$BYM2.C.prob15, carto$`Region`, function(x) sum(x>0.9))
Deaths.BYM2 <- tapply(carto[carto$BYM2.prob15>0.9,]$O, carto[carto$BYM2.prob15>0.9,]$`Region`, sum)
Deaths.BYM2.C <- tapply(carto[carto$BYM2.C.prob15>0.9,]$O, carto[carto$BYM2.C.prob15>0.9,]$`Region`, sum)
Overlapping.BYM2 <- table(carto[carto$BYM2.prob1>0.9 & carto$BYM2.nb.prob1>0.9,]$`Region`)
Overlapping.BYM2.C <- table(carto[carto$BYM2.C.prob1>0.9 & carto$BYM2.C.nb.prob1>0.9,]$`Region`)

Table8 <- cbind(Obs.SMR,High.BYM2,High.BYM2.C,Deaths.BYM2,Deaths.BYM2.C,Overlapping.BYM2,Overlapping.BYM2.C)
Table8[is.na(Table8)] <- 0
Table8 <- Table8[c(c("North East","North West","Yorkshire-Humberside","West Midlands","East Midlands","East","South East","London","South West")),]
Table8 <- rbind(Table8,c(1,apply(Table8,2,sum)[-1]))
Table8


################################################
## Table 9: Relative risk categories by model ##
################################################

## Compute posterior exceedence probabilities of being less than 0.67 (1/1.5) ##
carto$BYM2.prob.low <- unlist(lapply(BYM2.RR$marginals.fitted.values, function(x) inla.pmarginal(1/1.5, x)))
carto$BYM2.C.prob.low <- unlist(lapply(BYM2.C.RR$marginals.fitted.values, function(x) inla.pmarginal(1/1.5, x)))

Extreme.High <- c(sum(carto$BYM2.prob15>0.9),
                  round(100*mean(carto$BYM2.prob15>0.9),1),
                  sum(carto$BYM2.C.prob15>0.9),
                  round(100*mean(carto$BYM2.C.prob15>0.9),1))
High <- c(sum(carto$BYM2.prob1>0.9 & carto$BYM2.prob15<0.9),
          round(100*mean(carto$BYM2.prob1>0.9 & carto$BYM2.prob15<0.9),1),
          sum(carto$BYM2.C.prob1>0.9 & carto$BYM2.C.prob15<0.9),
          round(100*mean(carto$BYM2.C.prob1>0.9 & carto$BYM2.C.prob15<0.9),1))
Intermediate <- c(sum(carto$BYM2.prob1<0.9 & carto$BYM2.prob1>0.1),
                  round(100*mean(carto$BYM2.prob1<0.9 & carto$BYM2.prob1>0.1),1),
                  sum(carto$BYM2.C.prob1<0.9 & carto$BYM2.C.prob1>0.1),
                  round(100*mean(carto$BYM2.C.prob1<0.9 & carto$BYM2.C.prob1>0.1),1))
Low <- c(sum(carto$BYM2.prob1<0.1 & carto$BYM2.prob.low<0.9),
         round(100*mean(carto$BYM2.prob1<0.1 & carto$BYM2.prob.low<0.9),1),
         sum(carto$BYM2.C.prob1<0.1 & carto$BYM2.C.prob.low<0.9),
         round(100*mean(carto$BYM2.C.prob1<0.1 & carto$BYM2.C.prob.low<0.9),1))
Extreme.Low <- c(sum(carto$BYM2.prob.low>0.9),
                 round(100*mean(carto$BYM2.prob.low>0.9),1),
                 sum(carto$BYM2.C.prob.low>0.9),
                 round(100*mean(carto$BYM2.C.prob.low>0.9),1))

Table9 <- rbind(Extreme.High,High,Intermediate,Low,Extreme.Low)
colnames(Table9) <- c("BYM2 Model","% All MSOAs","BYM2+C Model","% All MSOAs")
Table9


###########################################################################
## Table A1: Model selection criteria for iCAR+C models fitted with INLA ##
##          considering different neighborhood orders (parameter l)      ##
###########################################################################
Table.A1 <- do.call(rbind,lapply(iCAR.C, DIC))
round(Table.A1,1)


###########################################################################
## Table A2: Model selection criteria for LCAR+C models fitted with INLA ##
##          considering different neighborhood orders (parameter l)      ##
###########################################################################
Table.A2 <- do.call(rbind,lapply(LCAR.C, DIC))
round(Table.A2,1)


###########################################################################
## Table A3: Model selection criteria for LCAR+C models fitted with INLA ##
##          considering different neighborhood orders (parameter l)      ##
###########################################################################
Table.A3 <- do.call(rbind,lapply(BYM.C, DIC))
round(Table.A3,1)
