library(MASS)
library(matrixcalc)
layered_sim2EigenCors<-function(layered_sim, ...)
{
    " Create correlation matrices between eigengenes nested in layered_sim object,
      with nested structure as one returned from LayeredSimulationByReplacement."
    Hier<- layered_sim$hierarchy
    hnames<- names(Hier)
        clist<-lapply(Hier, function(x) cor(x$MEs[,colnames(x$MEs)!='MEgrey'], ...))
    if (!is.null(hnames)){
        names(clist)<- hnames 
    }
    return(clist)
}

Cormat2Eigengenes<-function(ref_cor, n_samples) 
{
    if (!is.positive.definite(ref_cor)) stop("input correlation matrix should be positive definite")
    EIDF=as.data.frame(mvrnorm(n_samples, mu= rep(0, ncol(ref_cor)), Sigma=ref_cor))
    if (!is.null(colnames(ref_cor))) colnames(EIDF)<- colnames(ref_cor)
    return(EIDF)

}

layered_sim2
