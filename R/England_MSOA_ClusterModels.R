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


#####################################################################
## Apply DBSC algorithm over the residual of the non-spatial model ##
#####################################################################
if(!file.exists("DBSC_nbMatrices.Rdata")){
  
  X <- cbind(1,data.INLA[,c("ISOL","NH","HDD","AIRQ")])
  colnames(X)[1] <- "(Intercept)"
  head(X)
  
  carto$res <- log(carto$O/carto$E+0.0001)-as.matrix(X)%*%GLM$summary.fixed$`0.5quant`
  
  
  ## Compute l-order neighbourhood matrix ##
  l <- 8
  Wk <- vector("list",l)
  names(Wk) <- paste("W",seq(1,l),sep="")
  
  Wk[[1]] <- W
  for(i in seq(1,l-1)){
    Wk[[i+1]] <- Wk[[i]]%*%Wk[[1]]
    diag(Wk[[i+1]]) <- 0
    Wk[[i+1]][Wk[[i+1]]>0] <- 1
  }
  
  
  ## Run the DBSC algorithm (min.size=NULL) for each Wk ##
  Q.clust <- vector("list",l)
  names(Q.clust) <- paste("l",seq(l),sep="=")
  
  data.INLA.C <- rep(list(data.INLA),l)
  names(data.INLA.C) <- paste("l",seq(l),sep="=")
  
  for(i in seq(l)){
    cat("\nNeighborhood order: l=",i,"\n",sep="")
    
    cat("  + Running DBSC algorithm\n")
    carto.aux <- clustering_partition(carto=carto, ID.area="CODE", var="res", W=W, Wk=Wk[[i]],
                                      n.cluster=NULL, min.size=NULL, verbose=FALSE)
    
    cat("  + Computing neighbourhood graph of cluster-level random effect\n")
    k <- length(table(carto.aux$ID.group))
    carto.partition <- unionSpatialPolygons(as(carto.aux,"Spatial"),carto.aux$ID.group)
    carto.nb <- poly2nb(carto.partition)
    carto.W <- nb2mat(carto.nb, style="B")
    
    Q.clust[[i]] <- Diagonal(k,colSums(carto.W))-carto.W
    
    data.INLA.C[[i]]$ID.clust <- as.numeric(carto.aux$ID.group)
  }
  
  save(list=c("l","Wk","Q.clust","data.INLA.C"), file="DBSC_nbMatrices.Rdata")
  
}


##################################################################
## Define the hyperprior distributions and data for INLA models ##
##################################################################
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


#################################
## iCAR+C models (confounding) ##
#################################
load("DBSC_nbMatrices.Rdata")

iCAR.C <- vector("list",l)
names(iCAR.C) <- paste("l",seq(l),sep="=")

for(i in seq(l)){
  cat("Fitting INLA model: l=",i,"\n",sep="")
  
  formula <- O ~ ISOL + NH + HDD + AIRQ + 
    f(ID.area, model="besag", graph=Qs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
    f(ID.clust, model="besag", graph=Q.clust[[i]], constr=TRUE, hyper=list(prec=list(prior=sdunif)))

  Model <- inla(formula, family="poisson", data=data.INLA.C[[i]], E=E,
                control.predictor=list(compute=TRUE, cdf=c(log(1))),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                control.inla=list(strategy="simplified.laplace"))
  
  iCAR.C[[i]] <- Model
}

save(iCAR.C, file="England_MSOA_ClusterModels_iCAR.Rdata")


#################################
## LCAR+C models (confounding) ##
#################################
load("DBSC_nbMatrices.Rdata")

LCAR.C <- vector("list",l)
names(LCAR.C) <- paste("l",seq(l),sep="=")

for(i in seq(l)){
  cat("Fitting INLA model: l=",i,"\n",sep="")
  
  Q.clust.Leroux <- Diagonal(nrow(Q.clust[[i]]))-Q.clust[[i]]
  
  formula <- O ~ ISOL + NH + HDD + AIRQ + 
    f(ID.area, model="generic1", Cmatrix=Qs.Leroux, constr=TRUE,
      hyper=list(prec=list(prior=sdunif), beta=list(prior=lunif))) + 
    f(ID.clust, model="generic1", Cmatrix=Q.clust.Leroux, constr=TRUE,
      hyper=list(prec=list(prior=sdunif), beta=list(prior=lunif)))
  
  Model <- inla(formula, family="poisson", data=data.INLA.C[[i]], E=E,
                control.predictor=list(compute=TRUE, cdf=c(log(1))),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                control.inla=list(strategy="simplified.laplace"))
  
  LCAR.C[[i]] <- Model
}

save(LCAR.C, file="England_MSOA_ClusterModels_LCAR.Rdata")


################################
## BYM+C models (confounding) ##
################################
load("DBSC_nbMatrices.Rdata")

BYM.C <- vector("list",l)
names(BYM.C) <- paste("l",seq(l),sep="=")

for(i in seq(l)){
  cat("Fitting INLA model: l=",i,"\n",sep="")
  
  Q.clust.Leroux <- Diagonal(nrow(Q.clust[[i]]))-Q.clust[[i]]
  
  formula <- O ~ ISOL + NH + HDD + AIRQ + 
    f(ID.area, model="bym", graph=Qs, constr=TRUE,
      hyper=list(theta1=list(prior=sdunif), theta2=list(prior=sdunif))) + 
    f(ID.clust, model="bym", graph=Qs, constr=TRUE,
      hyper=list(theta1=list(prior=sdunif), theta2=list(prior=sdunif)))
  
  Model <- inla(formula, family="poisson", data=data.INLA.C[[i]], E=E,
                control.predictor=list(compute=TRUE, cdf=c(log(1))),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                control.inla=list(strategy="simplified.laplace"))
  
  BYM.C[[i]] <- Model
}

save(BYM.C, file="England_MSOA_ClusterModels_BYM.Rdata")


#################################
## BYM2+C models (confounding) ##
#################################
load("DBSC_nbMatrices.Rdata")

BYM2.C <- vector("list",l)
names(BYM2.C) <- paste("l",seq(l),sep="=")

for(i in seq(l)){
  cat("Fitting INLA model: l=",i,"\n",sep="")
  
  formula <- O ~ ISOL + NH + HDD + AIRQ + 
    f(ID.area, model="bym2", graph=Qs, constr=TRUE,
      hyper=list(prec=list(prior=sdunif), phi=list(prior=lunif, initial=0))) + 
    f(ID.clust, model="bym2", graph=Q.clust[[i]], constr=TRUE,
      hyper=list(prec=list(prior=sdunif), phi=list(prior=lunif, initial=0)))
  
  Model <- inla(formula, family="poisson", data=data.INLA.C[[i]], E=E,
                control.predictor=list(compute=TRUE, cdf=c(log(1))),
                control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE),
                control.inla=list(strategy="simplified.laplace"))
  
  BYM2.C[[i]] <- Model
}

save(BYM2.C, file="England_MSOA_ClusterModels_BYM2.Rdata")