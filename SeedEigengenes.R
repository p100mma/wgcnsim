library(MASS)

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
