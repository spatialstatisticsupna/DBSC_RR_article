rm(list=ls())
library(bigDM)
library(INLA)
library(maptools)
library(sf)
library(spdep)


#########################################
## Load the data and cartography files ##
#########################################
load("../data/England_MSOA.Rdata")
head(carto)

carto$HDD <- -carto$HDD ## Inverse measure

S <- nrow(carto)

Qs <- Diagonal(S,colSums(W))-W
Qs.Leroux <- Diagonal(S)-Qs

X <- as.matrix(st_set_geometry(carto, NULL)[,c("ISOL","NH","HDD","AIRQ")])
X <- apply(X,2,function(x) scale(x))
rownames(X) <- carto$CODE

data.INLA <- data.frame(O=carto$O, E=carto$E, X, Area=carto$CODE, ID.area=1:S)


###########################################################################################################
## Load previously fitted models. The BYM2 model fitted with INLA can be downloaded from:                ##
##  - https://emi-sstcdapp.unavarra.es/England_MSOA/INLA_models/England_MSOA_SpatialModels_BYM2.Rdata    ##
###########################################################################################################
load("DBSC_nbMatrices.Rdata")
load("England_MSOA_ClusterModels_BYM2.Rdata")

## Select l=1 neighborhood model ##
BYM2.C <- MODELS$`l=1`
summary(BYM2.C)

k <- length(table(data.INLA.C$`l=1`$ID.clust))
Q.clust <- Q.clust$`l=1`


#####################################################
## Fit the BYM2+C+RR model (restricted regression) ##
#####################################################
p <- ncol(X)
W <- Diagonal(S, BYM2.C$summary.fitted.values$mode*data.INLA$E)
W.sqrt <- Diagonal(S, sqrt(diag(W)))

XX <- cbind(rep(1,S),as.matrix(data.INLA[,c("ISOL","NH","HDD","AIRQ")]))
Pc <- Diagonal(S)-W.sqrt%*%XX%*%solve(t(XX)%*%W%*%XX)%*%t(XX)%*%W.sqrt
eigen.Pc <- eigen(Pc)
L <- eigen.Pc$vectors[,eigen.Pc$values>1e-12]
Z.area <- solve(W.sqrt)%*%L%*%t(L)%*%W.sqrt

aux <- fastDummies::dummy_cols(.data=data.INLA.C$`l=1`, select_columns="ID.clust")
aux <- as(aux[,seq(ncol(aux)-k+1,ncol(aux))],"Matrix")
Z.clust <- solve(W.sqrt)%*%L%*%t(L)%*%W.sqrt%*%aux

M0 <- solve(t(XX)%*%XX)%*%t(XX)
beta.lc = inla.make.lincombs(Predictor=M0, ID.area=-M0%*%Z.area, ID.clust=-M0%*%Z.clust)
names(beta.lc) <- paste("X",as.character(0:p),sep="")


## Define the hyperprior distributions and data for INLA models ##
sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"

lunif = "expression:
    a = 1;
    b = 1;
    beta = exp(theta)/(1+exp(theta));
    logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
    log_jacobian = log(beta*(1-beta));
    return(logdens+log_jacobian)"


## INLA model ##
formula <- O ~ ISOL + NH + HDD + AIRQ + 
  f(ID.area, model="bym2", graph=Qs, constr=TRUE,
    hyper=list(prec=list(prior=sdunif), phi=list(prior=lunif, initial=0))) + 
  f(ID.clust, model="bym2", graph=Q.clust, constr=TRUE,
    hyper=list(prec=list(prior=sdunif), phi=list(prior=lunif, initial=0)))

BYM2.C.RR <- inla(formula, family="poisson", data=data.INLA.C$`l=1`, E=E,
                  control.predictor=list(compute=TRUE, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                  lincomb=beta.lc,
                  control.inla=list(strategy="simplified.laplace"))

names <- rownames(BYM2.C.RR$summary.fixed)
BYM2.C.RR$summary.fixed <- BYM2.C.RR$summary.lincomb.derived[,-1]
rownames(BYM2.C.RR$summary.fixed) <- names
BYM2.C.RR$summary.lincomb.derived <- NULL

names <- names(BYM2.C.RR$marginals.lincomb.derived)
BYM2.C.RR$marginals.fixed <- BYM2.C.RR$marginals.lincomb.derived
names(BYM2.C.RR$marginals.fixed) <- names
BYM2.C.RR$marginals.lincomb.derived <- NULL

summary(BYM2.C.RR)

save(BYM2.C.RR, file="England_MSOA_ClusterModels_BYM2_RR.Rdata")
