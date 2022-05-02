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

CloneLayered<- function(ref_layered, ref_sim_argsList, eigen_method_cor='pearson', verbose=0,other_ref_simargsList=NULL)
{
    
    if (is.null(names(ref_sim_argsList))) stop("names of ref_sim_argsList must not be NULL")
    if (sum(names(ref_sim_argsList)!=names(ref_layered$hierarchy))>0) stop("names of ref_sim_argsList must be names of ref_layered$hierarchy")
    
    ClonedHier<- CloneLayeredHierarchy(ref_layered, eigen_method_cor)
    MergedResult=list(expr_data='foo',
                      hierarchy=ClonedHier,
                      GrayArea=ref_layered$GrayArea,
                      colors_final=ref_layered$colors_final)
    layer_counter=1
    for (layer in names(ref_sim_argsList))
    {
    other_args_=NULL
    if (!is.null(other_ref_simargsList)) other_args_=other_ref_simargsList[[layer]]
    args_<- ref_sim_argsList[[layer]]
    args_$eigengenes= ClonedHier[[layer]]$MEs
    if (layer_counter==1)
    { SIMRES=simulateDatExpr_fromInput(args_, verbose=verbose, other_named_args=other_args_)
      AdjSIMRES<- ReorderSimByReal(SIMRES$sim_result, list(color_labels=ClonedHier[[layer]]$color_labels) )
      MergedResult$expr_data<-AdjSIMRES$datExpr
    }
    else{
        temp_sim=simulateDatExpr_fromInput(args_, verbose=verbose, other_named_args=other_args_)
         AdjSIMRES<- ReorderSimByReal(temp_sim$sim_result, list(color_labels=ClonedHier[[layer]]$color_labels) )
        MergedResult$expr_data[,ClonedHier[[layer]]$color_labels!="grey"]=AdjSIMRES$datExpr[,ClonedHier[[layer]]$color_labels!="grey"]
        }
    layer_counter= layer_counter+1
    }
    return(MergedResult) 
}
