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

CloneLayeredSimEigengenes<- function(ref_layered,...)
{
    " From a layeredSim object like from LayeredSimulationByReplacement
      create a named list of eigengenes from N(0,1) but with cor like eigengenes
      in ref_layered (on each level) (no Grey Eigengene included) "
      n_samples<- nrow(ref_layered$expr_data)
      corlist<-layered_sim2EigenCors(ref_layered, ...)
      EGlist <-lapply(corlist, Cormat2Eigengenes, n_samples)
      if (!is.null(names(corlist))) names(EGlist)<- names(corlist)
      return(EGlist)
}

CloneLayeredHierarchy<- function(ref_layered, eigen_method_cor='pearson')
{
    Hier<- ref_layered$hierarchy
    if (is.null(names(Hier))) stop("names of ref_layered$hierarchy must not be NULL")
    EGlist<- CloneLayeredSimEigengenes(ref_layered, method=eigen_method_cor)
    COLORS<- lapply(Hier,function(x)  x$color_labels)
    names(COLORS)<- names(EGlist)
    clone<- lapply(seq_along(EGlist), function(i) list(color_labels=COLORS[[ names(EGlist)[[i]]  ]], MEs= EGlist[[ names( EGlist)[[i]]  ]]))
    names(clone)<-names(Hier)
    return(clone)
}

CloneLayered<- function(ref_layered, cormat, true_grey_frac=0.5, eigen_method_cor='pearson')
{
  neg_cor_prop <- cor2sim_neg_cor_prop(cormat)
  MEnames <- names(MEs)
  colornames<- sub(".*ME", "", MEnames)
  colornames<- colornames[ colornames!="grey"]     #skip grey module (independent genes)
  mod_props<- colors2props(color_labels, colornames, true_grey_frac)
  
final_args<- list(eigengenes= MEs[names(MEs)!="MEgrey"], nGenes=ncol(ref_layered$expr_data), 
                    modProportions=mod_props,
                    minCor=cor_par_list$min_cor, maxCor=cor_par_list$max_cor,
                    propNegativeCor=neg_cor_prop,
                    geneMeans=gmeans)
}
