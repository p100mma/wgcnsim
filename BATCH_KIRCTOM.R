source('initial.R')

cormat<- readRDS('KIRCmat.rds')

TOMmat<-TOMsimilarity(abs(cormat)^3)
saveRDS(TOMmat, 'KIRCTOM.rds')
