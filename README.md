# A clustering adjusted regression model to identify extreme risks in small areas
This repository contains the R code to fit with INLA the spatial models described in the work entitled _"Identifying extreme COVID-19 mortality risks in English small areas: a disease cluster approach"_ [(Adin et al., 2022)](10.21203/rs.3.rs-864393/v1). It also contains the necessary functions to reproduce all the figures and tables of the article.

In this work we consider several classical disease mapping models and an extension of the density-based spatial clustering (DBSC) algorithm proposed in Santafé et al. (2021) to an analysis of COVID-19 related mortality in English small areas during the first wave of the epidemic in the first half of 2020. 


## Table of contents

- [Data](#Data)
- [R code](#R-code)
- [References](#References)


# Data
COVID-19 mortality data and potential risks factors (covariates) across the 6791 middle super output areas (MSOAs) of England during the period March to June 2020. The data is publically available online without any form of restriction or copyright.

- The [**England_MSOA.Rdata**](https://github.com/spatialstatisticsupna/DBSC_RR_article/blob/master/data/England_MSOA.Rdata) file contains the following objects:
  	- **carto**: `sf` object containing the polygons of the MSOAs of England and 12 variables
  	- **W**: spatial adjacency matrix of the MSOAs of England
  

The COVID-19 deaths data is associated with the online article by the UK Office of National Statistics entitled _"Deaths involving COVID-19 by local area and socioeconomic deprivation: deaths occurring between 1 March and 31 July 2020"_ (Office of National Statistics (ONS), 2020). Data on ethnicity and nursing homes are from the UK Census, data on health deprivation are from a 2019 compendium of different types of small area deprivation (Ministry of Housing, Communities and Local Government (MHCLG), 2019), while data on air pollution are from the Access to Healthy Assets and Hazards small area indicators profile at https://www.cdrc.ac.uk/new-update-access-to-healthy-assets-and-hazards-ahah-data-resource/.

# R code


# Acknowledgements
This work has been supported by Projects MTM2017-82553-R (AEI/FEDER, UE) and PID2020-113125RB-I00/MCIN/AEI/10.13039/501100011033. 

![plot](https://github.com/spatialstatisticsupna/bigDM/blob/master/micin-aei.jpg)


# References
[Adin, A., Condong, P., Santafé, G., and Ugarte, M.D. (2022). Identifying extreme COVID-19 mortality risks in English small areas: a disease cluster approach. _Stochastic Environmental Research and Risk Assessment (under review)_. DOI: 10.21203/rs.3.rs-864393/v1](https://www.researchsquare.com/article/rs-864393/v1)

[Ministry of Housing, Communities and Local Government (MHCLG). (2019). _English Indices of Deprivation 2019_. MHCLG, London, UK](https://dera.ioe.ac.uk/34259/1/IoD2019_Technical_Report.pdf)

[Office of National Statistics (ONS). (2020). Deaths involving COVID-19 by local area and socioeconomic deprivation: deaths occurring between 1 March and 31 July 2020. _Statistical Bulletin_. ONS, London, UK](https://backup.ons.gov.uk/wp-content/uploads/sites/3/2020/08/Deaths-involving-COVID-19-by-local-area-and-socioeconomic-deprivation-deaths-occurring-between-1-March-and-31-.pdf)

[Santafé, G., Adin, A., Lee, D., and Ugarte, M.D. (2021). Dealing with risk discontinuities to estimate cancer mortality risks when the number of small areas is large. _Statistical Methods in Medical Research_, __30(1)__, 6-21.](https://doi.org/10.1177/0962280220946502) 
