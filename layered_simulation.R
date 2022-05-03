LayeredSimulationByReplacement <- function(baseSim, baseAssumedColors, baseAssumedMEs, 
                                           withinSims, withinAssumedColors, withinAssumedMEs,
                                           baseSimName="base",
                                           withinNames=paste0("within",seq_along(withinSims)))
{
  #Inputs: baseSim        <-   expression data simulated from real data using clustering it with...
  #     ...baseAssumedColors   and...
  #     ...baseAssumedMEs.     data is assumed to be ordered like real data up to ordering of modules.
  #     ...withinSims     <-   list of expression data simulated from real data using correspoding...
  #     ...withinAssumedColors ,list of membership vectors and...
  #     ...withinAssumedMEs    ,list of eigengenes data frames or alike. Again, ordered like real data.
  #        baseSimName and withinNames <- string and vector of strings, names of elements of hierarchy.
  #Result: list with elements:
  #     expr_data         <-   expression data after replacing baseSim columns with non-gray columns of...
  #                            ...withinSims.
  #     hierarchy         <-   named list, names(hierarchy)=c(baseSimName, withinNames), each element...
  #                            ...a list of 2: color_labels, MEs, which are corresponding original data...
  #                            ...for each layer.
  #     GrayArea          <-   logical conjunction of gray area indicators from all layers
  ALLgrey<- baseAssumedColors=="grey"
  RESexpr<- baseSim
  colors_final<- baseAssumedColors
  withinHierarchy<- vector(mode="list", length=length(withinAssumedColors))
  for (i in seq_along(withinAssumedColors))
      { ALLgrey<- (ALLgrey & (withinAssumedColors[[i]]=="grey") ) 
        RESexpr[,withinAssumedColors[[i]]!="grey"]= withinSims[[i]][,withinAssumedColors[[i]]!="grey"]
        #^ take non-grey modules from withinSims, put in resulting simulation at their placement
        colors_final[withinAssumedColors[[i]]!="grey"]=paste0(withinAssumedColors[[i]][withinAssumedColors[[i]]!="grey"],i)
        withinHierarchy[[i]]<- list(color_labels= withinAssumedColors[[i]], MEs= withinAssumedMEs[[i]])
      }
  hierarchy=c(list(list(color_labels=baseAssumedColors, MEs= baseAssumedMEs)), withinHierarchy)
  names(hierarchy)= c(baseSimName, withinNames)
  return( list( 
                expr_data=RESexpr,
                hierarchy=hierarchy,
                GrayArea=ALLgrey,
                colors_final=colors_final
              )
        )
}

SmallGrayModules<- function(GrayIndicator, sim_expr_data, Other_full_clustering, neg_cor_prop,
                            MAX_COR=1){
  #1. Find clusters in Other_full_clustering bigger than 1 and in Gray area
  #2. For every such cluster in Other_full_clustering get its IDX
  #3. For every such cluster simulate module from random noise variable
  #4. Put every simulated module in the corresponding IDX in sim_expr_data
  counts<-table(Other_full_clustering[GrayIndicator])
  nonsingles_labels<-names(counts>1)
  counter=1
  print(paste0("there are ",length(nonsingles_labels),"clusters to fill."))
  for (label in nonsingles_labels){
    print(paste0("During cluster ",counter))
    LabeledAndGray<- (Other_full_clustering==label)&(GrayIndicator)
    modSize<- sum(LabeledAndGray)
    minCor=runif(1,0,MAX_COR)
    maxCor=runif(1, minCor,MAX_COR)
    ME<- rnorm(nrow(sim_expr_data))
    MOD<-simulateModule(ME, modSize, 
      nNearGenes = 0, 
      minCor =minCor , maxCor = maxCor, corPower = 1, 
      signed = FALSE, propNegativeCor = neg_cor_prop, 
      geneMeans = NULL,
      verbose = 0, indent = 0)
    sim_expr_data[,LabeledAndGray] = MOD
    counter=counter+1
  }
  return(sim_expr_data)
}

RandomGreyModules<- function(GrayIndicator, sim_expr_data, neg_cor_prop,
                             MAX_COR=0.75, MIN_COR=0, MaxSubmoduleSize=5, corpower=1,
                             sizeProbs=rep( 1/MaxSubmoduleSize,MaxSubmoduleSize)){
  ColorIDX<- which(GrayIndicator==TRUE)
  counter=1
  modSize<-sample(1:MaxSubmoduleSize, 1, prob=sizeProbs)
  while(counter+modSize<sum(GrayIndicator)){
    minCor=runif(1,MIN_COR,MAX_COR)
    maxCor=runif(1, minCor,MAX_COR)
    ME<- rnorm(nrow(sim_expr_data))
    MOD<-simulateModule(ME, modSize,
                        nNearGenes=0,
                        minCor=minCor,maxCor=maxCor, corPower=corpower,
                        signed=FALSE, propNegativeCor=neg_cor_prop,
                        geneMeans=NULL,
                        verbose=0, indent=0)
    sim_expr_data[, ColorIDX[counter:(counter+modSize-1)]] = MOD
    modSize<-sample(1:MaxSubmoduleSize, 1)
    counter= counter + modSize
  }
  return(sim_expr_data)
}


TestRandomGreyModules<- function(GrayIndicator, sim_expr_data, neg_cor_prop,
                                 MAX_CORs, MIN_CORs, MaxSubmoduleSizes,
                                 sizeProbsList=NULL, method_cor="spearman", 
                                 powerEstimate=5) {
  #Test for mean degree in grey area as a fun of MAX_COR, MIN_COR, MaxSubmoduleSize.
  mean_degs<- vector(length=length(MaxSubmoduleSizes))
  for (i in seq_along(MaxSubmoduleSizes)){
    print(paste0(i,"/",length(MaxSubmoduleSizes)," parameters."))
    print(paste0("MAX: ",MAX_CORs[[i]],"MIN: ",MIN_CORs[[i]],"MaxModSize: ",MaxSubmoduleSizes[[i]]))
    if(is.null(sizeProbsList)){
      PROBS= rep(1/MaxSubmoduleSizes[[i]], MaxSubmoduleSizes[[i]])
    } else { PROBS=sizeProbsList[[i]]  }
    tempExpr <- RandomGreyModules(GrayIndicator, sim_expr_data, neg_cor_prop,
                                  MAX_CORs[[i]], MIN_CORs[[i]], MaxSubmoduleSizes[[i]],
                                  PROBS)
    tempAdj<- abs(cor(tempExpr, method= method_cor))^powerEstimate
    degree_distrib<- intramodularConnectivity(tempAdj, rep("blue", ncol(tempExpr)))
    degree_distrib<- degree_distrib$kTotal
    mean_degs[[i]]= mean(degree_distrib[GrayIndicator])
    print(paste0("mean_deg: ",mean_degs[[i]]))
  }
  return(mean_degs)
}

SmallSubModules<- function(base_color_labels, sim_expr_data, Other_full_clustering, neg_cor_prop,
                           MAX_COR=0.625, MIN_COR=0.375, MaxSubmoduleSize=40){
  all_colors<- unique(base_color_labels)
  for (color_label in all_colors){if(color_label!="grey"){
    ColorIndicator<- base_color_labels==color_label
    counts<-table(Other_full_clustering[ColorIndicator])
    OKsizes_labels<- names(counts[(counts>1)&(counts<MaxSubmoduleSize)])
    print(paste0("there are ",length(OKsizes_labels),"clusters to fill in module ",color_label,"."))
    for (label in OKsizes_labels){
      LabeledAndCurrentColor<- (Other_full_clustering==label)&(ColorIndicator)
      modSize<- sum(LabeledAndCurrentColor)
      minCor=runif(1,MIN_COR,MAX_COR)
      maxCor=runif(1, minCor,MAX_COR)
      ME<- colors2MEs(rep("blue", modSize), as.data.frame(sim_expr_data[, LabeledAndCurrentColor]))
      ME<- ME$eigengenes[,1]
      MOD<-simulateModule(ME, modSize,
                          nNearGenes=0,
                          minCor=minCor,maxCor=maxCor, corPower=1,
                          signed=FALSE, propNegativeCor=neg_cor_prop,
                          geneMeans=NULL,
                          verbose=0, indent=0)
      sim_expr_data[,LabeledAndCurrentColor] = MOD
    }
    
  }}
  return(sim_expr_data)
}

RandomSmallSubModules<- function(base_color_labels, sim_expr_data, neg_cor_prop,
                                 MAX_COR=0.625, MIN_COR=0.375, MaxSubmoduleSize=40,
                                 sizeProbs=rep( 1/MaxSubmoduleSize,MaxSubmoduleSize)){
  all_colors<- unique(base_color_labels)
  for (color_label in all_colors){if(color_label!="grey"){
    print(paste0("Filling ",color_label))
    ColorIndicator<- base_color_labels==color_label
    ColorIDX<- which(base_color_labels==color_label)
    counter=1
    modSize<-sample(1:MaxSubmoduleSize, 1, prob=sizeProbs)
    while(counter+modSize<sum(ColorIndicator)){
      minCor=runif(1,MIN_COR,MAX_COR)
      maxCor=runif(1, minCor,MAX_COR)
      ME<- colors2MEs(rep("blue", modSize), as.data.frame(sim_expr_data[, ColorIDX[counter:(counter+modSize-1)]]))
      ME<- ME$eigengenes[,1]
      MOD<-simulateModule(ME, modSize,
                          nNearGenes=0,
                          minCor=minCor,maxCor=maxCor, corPower=1,
                          signed=FALSE, propNegativeCor=neg_cor_prop,
                          geneMeans=NULL,
                          verbose=0, indent=0)
      sim_expr_data[, ColorIDX[counter:(counter+modSize-1)]] = MOD
      modSize<-sample(1:MaxSubmoduleSize, 1)
      counter= counter + modSize
                                              }
    }}
    return(sim_expr_data)
  }


Layered2Submodules <- function(layered_simulation,  neg_cor_prop, base_sim_name='base',
                               GrMAX_COR=0.75, GrMIN_COR=0, GrMaxSubmoduleSize=5, Grcorpower=1,
                               GrsizeProbs=rep( 1/GrMaxSubmoduleSize,GrMaxSubmoduleSize),
                                ModMAX_COR=0.625, ModMIN_COR=0.375, ModMaxSubmoduleSize=40,
                                ModsizeProbs=rep( 1/ModMaxSubmoduleSize,ModMaxSubmoduleSize)){

layered_simulation$grayed<-RandomGreyModules(layered_simulation$GrayArea, layered_simulation$expr_data,
                             neg_cor_prop,
                             GrMAX_COR,GrMIN_COR, GrMaxSubmoduleSize, Grcorpower,
                             GrsizeProbs) 
layered_simulation$submods<-RandomSmallSubModules(layered_simulation$hierarchy[[base_sim_name]]$color_labels,
                                                 layered_simulation$grayed, neg_cor_prop,
                                                 ModMAX_COR, ModMIN_COR, ModMaxSubmoduleSize,
                                                 ModsizeProbs)
return(layered_simulation)
    




}
# depends on the igraph library sizeProbs=rep( 1/MaxSubmoduleSize,MaxSubmoduleSize)
#
# TestRandomSmallSubModules <- function(base_color_labels, sim_expr_data, neg_cor_prop,
#                                       MaxSubmoduleSizes, sizeProbsList,MAX_CORs, MIN_CORs,
#                                       powerEstimate=5, method_cor="spearman", referenceSize=500
#                                       ){
#   ARIs=vector(length=length(MaxSubmoduleSizes))
#   VIs=vector(length=length(MaxSubmoduleSizes))
#   NMIs=vector(length=length(MaxSubmoduleSizes))
#   for (i in seq_along(MaxSubmoduleSizes)){
#     print(i)
#     TEMPDATA<- RandomSmallSubModules(base_color_labels, sim_expr_data, neg_cor_prop, MAX_CORs[[i]], MIN_CORs[[i]],
#                                      MaxSubmoduleSizes[[i]],sizeProbsList[[i]])
#     TEMPADJ<- abs(cor(TEMPDATA, method=method_cor))^powerEstimate
#     TEMPCLRES<- ClusteringResults.fromAdjacency(TEMPADJ, TEMPDATA, save_final=FALSE, minModSize=referenceSize,
#                                                 dc_doPAM=FALSE, calculateMEs=FALSE)
#     ARIs[[i]]<-compare(base_color_labels,TEMPCLRES$color_labels,"adjusted.rand")
#     VIs[[i]]<-compare(base_color_labels,TEMPCLRES$color_labels,"vi")
#     NMIs[[i]]<-compare(base_color_labels,TEMPCLRES$color_labels,"nmi")
#     dev_closure( jpeg, 
#                  {plotDendroAndColors(TEMPCLRES$geneTree, 
#                                       as.data.frame(list(TEMPCLRES$color_labels,base_color_labels)),
#                                       c(paste0("params",i),paste0("HS",referenceSize)),
#                                                     main=paste0("params",i),
#                                       hang=0.03, dendroLabels=FALSE
#                                       )}, 
#                  paste0("params",i,".jpg"), width=1000, height=600)
#   }
#   return(list(VIs=VIs,ARIs=ARIs,NMIs=NMIs,MaxSubmoduleSizes=MaxSubmoduleSizes,sizeProbsList=sizeProbsList,
#               MAX_CORs=MAX_CORs, MIN_CORs=MIN_CORs))
#   }
# 
# 
# TestRandomSmallSubModulesAndNetworkConcepts<- function(base_color_labels, sim_expr_data, neg_cor_prop,
#                                                        MaxSubmoduleSizes, sizeProbsList,MAX_CORs, MIN_CORs,
#                                                       powerEstimate=5, method_cor="spearman",
#                                                       referenceSize=500){
#   ARIs=vector(length=length(MaxSubmoduleSizes))
#   VIs=vector(length=length(MaxSubmoduleSizes))
#   NMIs=vector(length=length(MaxSubmoduleSizes))
#   Concepts=vector(mode="list",length=length(MaxSubmoduleSizes))
#   for (i in seq_along(MaxSubmoduleSizes)){
#     print(i)
#     TEMPDATA<- RandomSmallSubModules(base_color_labels, sim_expr_data, neg_cor_prop, MAX_CORs[[i]], MIN_CORs[[i]],
#                                      MaxSubmoduleSizes[[i]],sizeProbsList[[i]])
#     TEMPADJ<- abs(cor(TEMPDATA, method=method_cor))^powerEstimate
#     TEMPCLRES<- ClusteringResults.fromAdjacency(TEMPADJ, TEMPDATA, save_final=FALSE, minModSize=referenceSize,
#                                                 dc_doPAM=FALSE, calculateMEs=FALSE)
#     ARIs[[i]]<-compare(base_color_labels,TEMPCLRES$color_labels,"adjusted.rand")
#     VIs[[i]]<-compare(base_color_labels,TEMPCLRES$color_labels,"vi")
#     NMIs[[i]]<-compare(base_color_labels,TEMPCLRES$color_labels,"nmi")
#     dev_closure( jpeg, 
#                  {plotDendroAndColors(TEMPCLRES$geneTree, 
#                                       as.data.frame(list(TEMPCLRES$color_labels,base_color_labels)),
#                                       c(paste0("params_",i),paste0("HS",referenceSize)),
#                                       main=paste0("params_",i),
#                                       hang=0.03, dendroLabels=FALSE
#                  )}, 
#                  paste0("params_",i,".jpg"), width=1000, height=600)
#     print("Concepts...")
#     TEMPCEPTS<-fundamentalNetworkConcepts(TEMPADJ)
#     Concepts[[i]]<- TEMPCEPTS
#   }
#   return(list(VIs=VIs,ARIs=ARIs,NMIs=NMIs,MaxSubmoduleSizes=MaxSubmoduleSizes,sizeProbsList=sizeProbsList,
#               MAX_CORs=MAX_CORs, MIN_CORs=MIN_CORs, Concepts= Concepts))
# }
# 
# 
