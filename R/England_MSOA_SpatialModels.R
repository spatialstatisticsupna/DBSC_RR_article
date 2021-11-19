rm(list=ls())
library(bigDM)
library(INLA)
library(sf)


#########################################
## Load the data and cartography files ##
#########################################
load("../data/England_MSOA.Rdata")
head(carto)

carto$HDD <- -carto$HDD ## Inverse measure

S <- nrow(carto)

Qs <- Diagonal(S,colSums(W))-W
Qs.Leroux <- Diagonal(S)-Qs


#############################
## GLM (non-spatial model) ##
#############################
X <- as.matrix(st_set_geometry(carto, NULL)[,c("ISOL","NH","HDD","AIRQ")])
X <- apply(X,2,function(x) scale(x))
rownames(X) <- carto$CODE

data.INLA <- data.frame(O=carto$O, E=carto$E, X, Area=carto$CODE, ID.area=1:S)

f.GLM <- O ~ ISOL + NH + HDD + AIRQ

GLM <- inla(f.GLM, family="poisson", data=data.INLA, E=E,
            control.predictor=list(compute=TRUE, cdf=c(log(1))),
            control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
            control.inla=list(strategy="simplified.laplace"))

summary(GLM)


############################################################
## iCAR (confounding) and iCAR+RR (restricted regression) ##
############################################################
iCAR <- CAR_INLA(carto=carto, ID.area="CODE", O="O", E="E",
                 X=c("ISOL","NH","HDD","AIRQ"), confounding=NULL,
                 prior="intrinsic", W=W, model="global", strategy="simplified.laplace")

summary(iCAR)


iCAR.RR <- CAR_INLA(carto=carto, ID.area="CODE", O="O", E="E",
                 X=c("ISOL","NH","HDD","AIRQ"), confounding="restricted",
                 prior="intrinsic", W=W, model="global", strategy="simplified.laplace")

summary(iCAR.RR)

save(list=c("GLM","iCAR","iCAR.RR"), file="England_MSOA_SpatialModels_iCAR.Rdata")


############################################################
## LCAR (confounding) and LCAR+RR (restricted regression) ##
############################################################
LCAR <- CAR_INLA(carto=carto, ID.area="CODE", O="O", E="E",
                 X=c("ISOL","NH","HDD","AIRQ"), confounding=NULL,
                 prior="Leroux", W=W, model="global", strategy="simplified.laplace")

summary(LCAR)


LCAR.RR <- CAR_INLA(carto=carto, ID.area="CODE", O="O", E="E",
                    X=c("ISOL","NH","HDD","AIRQ"), confounding="restricted",
                    prior="Leroux", W=W, model="global", strategy="simplified.laplace")

summary(LCAR.RR)

save(list=c("GLM","LCAR","LCAR.RR"), file="England_MSOA_SpatialModels_LCAR.Rdata")


##########################################################
## BYM (confounding) and BYM+RR (restricted regression) ##
##########################################################
BYM <- CAR_INLA(carto=carto, ID.area="CODE", O="O", E="E",
                X=c("ISOL","NH","HDD","AIRQ"), confounding=NULL,
                prior="BYM", W=W, model="global", strategy="simplified.laplace")

summary(BYM)


BYM.RR <- CAR_INLA(carto=carto, ID.area="CODE", O="O", E="E",
                   X=c("ISOL","NH","HDD","AIRQ"), confounding="restricted",
                   prior="BYM", W=W, model="global", strategy="simplified.laplace")

summary(BYM.RR)

save(list=c("GLM","BYM","BYM.RR"), file="England_MSOA_SpatialModels_BYM.Rdata")


############################################################
## BYM2 (confounding) and BYM2+RR (restricted regression) ##
############################################################
BYM2 <- CAR_INLA(carto=carto, ID.area="CODE", O="O", E="E",
                 X=c("ISOL","NH","HDD","AIRQ"), confounding=NULL,
                 prior="BYM2", W=W, model="global", strategy="simplified.laplace")

summary(BYM2)


BYM2.RR <- CAR_INLA(carto=carto, ID.area="CODE", O="O", E="E",
                    X=c("ISOL","NH","HDD","AIRQ"), confounding="restricted",
                    prior="BYM2", W=W, model="global", strategy="simplified.laplace")

summary(BYM2.RR)

save(list=c("GLM","BYM2","BYM2.RR"), file="England_MSOA_SpatialModels_BYM2.Rdata")
