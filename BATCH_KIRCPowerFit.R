library(WGCNA)
source('inital.R')

cormat<- readRDS('KIRCmat.rds')

Result<-powerEstimation(cormat)
saveRDS(Result, 'KIRCpowerFit.rds')
