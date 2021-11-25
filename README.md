# A clustering adjusted regression model to identify extreme risks in small areas
This repository contains the R code to fit with INLA the spatial models described in the work entitled _"Identifying extreme COVID-19 mortality risks in English small areas: a disease cluster approach"_ [(Adin et al., 2022)](https://doi.org/10.21203/rs.3.rs-864393/v1). It also contains the necessary functions to reproduce all the figures and tables of the article.

In this work we consider several classical disease mapping models and an extension of the density-based spatial clustering (DBSC) algorithm proposed in Santafé et al. (2021) to an analysis of COVID-19 related mortality in English small areas during the first wave of the epidemic in the first half of 2020. 


## Table of contents

- [Data](#Data)
- [R code](#R-code)
- [References](#References)


# Data
COVID-19 mortality data and potential risks factors (covariates) across the 6791 middle super output areas (MSOAs) of England during the period March to June 2020. The data are publicly available online without any form of restriction or copyright.

The [**England_MSOA.Rdata**](https://github.com/spatialstatisticsupna/DBSC_RR_article/blob/master/data/England_MSOA.Rdata) file contains the following objects:
  - **W**: spatial adjacency matrix of the MSOAs of England
  - **carto**: `sf` object containing the polygons of the MSOAs of England and 12 variables
    - **_CODE_**: character vector of geographic identifiers
    - **_NAME_**: character vector of MSOAs name
    - **_Region_**: classification of MSOAs by region (factor with 9 levels)
    - **_Urban-Rural_**: classification of MSOAs according to urbal-rural category (factor with 8 levels)
    - **_O_**: observed number of COVID-19 related deaths
    - **_E_**: expected number of COVID-19 related deaths
    - **_SMR_**: standardized mortality ratio
    - **_ISOL_**: Lieberson isolation index
    - **_NH_**: nursing home location
    - **_HDD_**: health deprivation and disability index
    - **_AIRQ_**: measure of poor air quality
    - **_geometry_**: sfc_GEOMETRY
  

The COVID-19 deaths data is associated with the online article by the UK Office of National Statistics entitled _"Deaths involving COVID-19 by local area and socioeconomic deprivation: deaths occurring between 1 March and 31 July 2020"_ (Office of National Statistics (ONS), 2020). Data on ethnicity and nursing homes are from the UK Census, data on health deprivation are from a 2019 compendium of different types of small area deprivation (Ministry of Housing, Communities and Local Government (MHCLG), 2019), while data on air pollution are from the Access to Healthy Assets and Hazards small area indicators profile at https://www.cdrc.ac.uk/new-update-access-to-healthy-assets-and-hazards-ahah-data-resource/.


# R code
R code to fit the spatial models with INLA (http://www.r-inla.org/) considered in the present paper, and code to reproduce all the figures and tables. All the R files are written by the authors of the paper.

- [**England_MSOA_SpatialModels.R**](https://github.com/spatialstatisticsupna/DBSC_RR_article/blob/master/R/England_MSOA_SpatialModels.R)

  R code to fit spatial models using the [**bigDM**](https://github.com/spatialstatisticsupna/bigDM) package that incorporates area-level random effects (see Equation (1) of Adin et al., 2021), denoted as iCAR, LCAR, BYM or BYM2 depending on the conditional autoregressive (CAR) prior distribution considered for the spatially structured random effect. The restricted regression approach (Reich et al., 2006) has been adopted to deal with spatial confounding between fixed and random effects. For those models, the "+RR" suffix has been added to the notation.

- [**England_MSOA_ClusterModels.R**](https://github.com/spatialstatisticsupna/DBSC_RR_article/blob/master/R/England_MSOA_ClusterModels.R)

  R code to fit spatial models that incorporate both area-level and cluster-level random effects (see Equation (3) of Adin et al., 2022), denoted as iCAR+C, LCAR+C, BYM+C or BYM2+C depending on the CAR prior distribution adopted for both spatial random effects. Our modelling approach consists of a two-stage procedure. First, the DBSC algorithm (implemented in the [**bigDM**](https://github.com/spatialstatisticsupna/bigDM) package) is applied over the residuals of the non-spatial model (GLM). Then, a model with area-level and cluster-level random effects is fitted.
  
  *_NOTE: Adjacency matrices (`Wk`), neighborhood structure matrices (`Q.clust`) and final INLA datasets (`data.INLA.C`) that are necessary to fit these models have been previously precomputed and stored at [**DBSC_nbMatrices.Rdata**](https://github.com/spatialstatisticsupna/DBSC_RR_article/blob/master/R/DBSC_nbMatrices.Rdata) file to save computation time._

- [**England_MSOA_ClusterModels_BYM2_RR.R**](https://github.com/spatialstatisticsupna/DBSC_RR_article/blob/master/R/England_MSOA_ClusterModels_BYM2_RR.R)
  
  R code to fit the clustering adjusted regression model with the BYM2 prior distribution and restricted regression to deal with spatial confounding issues.

- [**Figures_and_Tables.R**](https://github.com/spatialstatisticsupna/DBSC_RR_article/blob/master/R/Figures_and_Tables.R)
  
  This R script contains the necessary functions to reproduce all the figures and tables of the present paper. The fitted models with INLA can be downloaded from https://emi-sstcdapp.unavarra.es/England_MSOA/INLA_models/.
  
  
# Acknowledgements
This work has been supported by Projects MTM2017-82553-R (AEI/FEDER, UE) and PID2020-113125RB-I00/MCIN/AEI/10.13039/501100011033. 

![plot](https://github.com/spatialstatisticsupna/bigDM/blob/master/micin-aei.jpg)


# References
[Adin, A., Congdon, P., Santafé, G., and Ugarte, M.D. (2022). Identifying extreme COVID-19 mortality risks in English small areas: a disease cluster approach. _Stochastic Environmental Research and Risk Assessment (under review)_. DOI: 10.21203/rs.3.rs-864393/v1](https://www.researchsquare.com/article/rs-864393/v1)

[Ministry of Housing, Communities and Local Government (MHCLG). (2019). _English Indices of Deprivation 2019_. MHCLG, London, UK](https://dera.ioe.ac.uk/34259/1/IoD2019_Technical_Report.pdf)

[Office of National Statistics (ONS). (2020). Deaths involving COVID-19 by local area and socioeconomic deprivation: deaths occurring between 1 March and 31 July 2020. _Statistical Bulletin_. ONS, London, UK](https://backup.ons.gov.uk/wp-content/uploads/sites/3/2020/08/Deaths-involving-COVID-19-by-local-area-and-socioeconomic-deprivation-deaths-occurring-between-1-March-and-31-.pdf)

[Reich, B.J., Hodges, J.S., Zadnik, V. (2006). Effects of residual smoothing on the posterior of the fixed effects in disease-mapping models. _Biometrics_, **62(4)**:1197-1206.](https://doi.org/10.1111/j.1541-0420.2006.00617.x)

[Santafé, G., Adin, A., Lee, D., and Ugarte, M.D. (2021). Dealing with risk discontinuities to estimate cancer mortality risks when the number of small areas is large. _Statistical Methods in Medical Research_, __30(1)__, 6-21.](https://doi.org/10.1177/0962280220946502) 
