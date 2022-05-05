
#TO DO: add option to evaluate min max TOM on Gray Submodules

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
        #print(dim(withinAssumedColors[[i]]))
        #print(dim(RESexpr))
        #print(dim(withinSims[[i]]))
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

SimSpecs2Hierarchy<- function( BaseSimSpec, WithinSimSpecs,
                             datasets_path,
                     base_expr_data_path=NULL, base_expr_data=NULL,
                     base_has_decision=FALSE, base_expr_RData=TRUE, method_cor='spearman', verbose=0){
   Hierarchy<- vector(mode='list', length=(length(WithinSimSpecs)+1))
   names(Hierarchy) <- c('base', paste0('within', seq_along(WithinSimSpecs)))
tempClres<-DoClusterFromFilenameArgs(datasets_path, BaseSimSpec$base_dataset_name, BaseSimSpec$base_network_name,BaseSimSpec$base_clustering_name, method_cor,
                             base_expr_data, base_expr_data_path, base_has_decision, base_expr_RData, 
calculateMEs=TRUE)
   Hierarchy$base<- list(color_labels=tempClres$color_labels, MEs=tempClres$ME_data$eigengenes, base_sim_specs=BaseSimSpec)
   for (i in seq_along(WithinSimSpecs)){    
    tempClres<-DoClusterFromFilenameArgs(datasets_path, WithinSimSpecs[[i]]$base_dataset_name, WithinSimSpecs[[i]]$base_network_name,WithinSimSpecs[[i]]$base_clustering_name, method_cor,
                             base_expr_data, base_expr_data_path, base_has_decision, base_expr_RData, 
calculateMEs=TRUE)
        Hierarchy[[paste0('within',i)]]<- list(color_labels=tempClres$color_labels, MEs=tempClres$ME_data$eigengenes, base_sim_specs=WithinSimSpecs[[i]])
}
return(Hierarchy)
}

SpecsFiles2Hierarchy<- function(BaseSimName, WithinSimNamesList, datasets_path, dataset_name,
                           base_expr_data_path=NULL, base_expr_data=NULL, #one of those must not be null
                           base_has_decision=FALSE, base_expr_RData=TRUE,
                           method_cor='spearman',verbose=0) 
{
    BaseSimSpec<- ReadSimSpecsFile(BaseSimName, datasets_path, dataset_name)
    WithinSimSpecs<- lapply(WithinSimNamesList, ReadSimSpecsFile, datasets_path, dataset_name)

Hierarchy<- SimSpecs2Hierarchy( BaseSimSpec, WithinSimSpecs,datasets_path,
                                base_expr_data_path, base_expr_data,
                     base_has_decision, base_expr_RData, method_cor, verbose)
return(Hierarchy)
}

SaveHierarchy<- function(Hierarchy, hierarchy_name, datasets_path) {
    dataset_name<- Hierarchy$base$base_sim_specs$base_dataset_name
    network_name<- Hierarchy$base$base_sim_specs$base_network_name
    prefix_path=paste0(datasets_path,'/',dataset_name,'/simulations/hierarchical/')
    dir.create(paste0(prefix_path, hierarchy_name), showWarnings=FALSE, recursive=TRUE)
    saveRDS(Hierarchy,paste0(prefix_path,hierarchy_name,'/hierarchy.rds'))
        }
ReadHierarchy<- function(hierarchy_name, datasets_path, dataset_name) {
    prefix_path=paste0(datasets_path,'/',dataset_name,'/simulations/hierarchical/')
    loadedH<-readRDS(paste0(prefix_path,hierarchy_name,'/hierarchy.rds'))
    return(loadedH)}
LayeredFromHierarchy<-function(Hierarchy, datasets_path,
                           base_expr_data_path=NULL, base_expr_data=NULL, #one of those must not be null
                           base_has_decision=FALSE, base_expr_RData=TRUE,
                           method_cor='spearman',verbose=0
                                )
{       
   baseAssumedColors<- Hierarchy$base$color_labels
   baseAssumedMEs<- Hierarchy$base$MEs
   baseSim<-Specs2Sim(Hierarchy$base$base_sim_spec,datasets_path,
                           base_expr_data_path, base_expr_data, #one of those must not be null
                           base_has_decision, base_expr_RData,
                           method_cor,verbose)$datExpr  
   withinAssumedColors<- vector(mode='list', length= length(Hierarchy)-1)
   withinAssumedMEs<- vector(mode='list', length= length(Hierarchy)-1)
   withinSims<- vector(mode='list', length= length(Hierarchy)-1)
    if (verbose>0) print("Done base")
   for (i in 1:(length(Hierarchy)-1)){
    withinAssumedColors[[i]]<- Hierarchy[[paste0('within',i)]]$color_labels
    withinAssumedMEs[[i]]<- Hierarchy[[paste0('within',i)]]$MEs
    withinSims[[i]]<- Specs2Sim(Hierarchy[[paste0('within',i)]]$base_sim_spec,datasets_path,
                           base_expr_data_path, base_expr_data, #one of those must not be null
                           base_has_decision, base_expr_RData,
                           method_cor,verbose)$datExpr
     
    if (verbose>0) print(paste0("Done within",i))
    }
   layered_sim<-LayeredSimulationByReplacement(baseSim, baseAssumedColors, baseAssumedMEs, 
                                           withinSims, withinAssumedColors, withinAssumedMEs)
return(layered_sim)
}
LayeredFromHierPath<-function(hierarchy_name, datasets_path,dataset_name,
                           base_expr_data_path=NULL, base_expr_data=NULL, #one of those must not be null
                           base_has_decision=FALSE, base_expr_RData=TRUE,
                           method_cor='spearman',verbose=0){

    prefix_path=paste0(datasets_path,'/',dataset_name,'/simulations/hierarchical/')
    Hier<-readRDS(paste0(prefix_path,hierarchy_name,'/hierarchy.rds'))
   return(LayeredFromHierarchy(Hier, datasets_path,
                           base_expr_data_path, base_expr_data, #one of those must not be null
                           base_has_decision, base_expr_RData,
                           method_cor,verbose)) 
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

RandomGreyModulesFromHier<- function(Hierarchy, datasets_path,MAX_COR,MIN_COR,MaxSubmoduleSize, neg_cor_prop,
                           corpower=1, sizeProbs=rep( 1/MaxSubmoduleSize,MaxSubmoduleSize),
                           base_expr_data_path=NULL, base_expr_data=NULL, #one of those must not be null
                           base_has_decision=FALSE, base_expr_RData=TRUE,
                           method_cor='spearman',verbose=0){
            
    LS<-LayeredFromHierarchy(Hierarchy, datasets_path,
                        base_expr_data_path, base_expr_data, #one of those must not be null
                        base_has_decision, base_expr_RData,
                        method_cor,verbose)
   LS$RGM<- RandomGreyModules(LS$GrayArea, LS$expr_data, neg_cor_prop,
                             MAX_COR, MIN_COR, MaxSubmoduleSize, corpower,
                             sizeProbs) 
    return(LS)
}

SaveRandomGreyModulesFromHierPath<- function(hierarchy_name, datasets_path,dataset_name,neg_cor_prop,
                                            GraySubmodulesName, MAX_COR,MIN_COR,MaxSubmoduleSize,
                                            corpower=1, sizeProbs=rep( 1/MaxSubmoduleSize,MaxSubmoduleSize),
                                            base_expr_data_path=NULL, base_expr_data=NULL, #one of those must not be null
                                            base_has_decision=FALSE, base_expr_RData=TRUE,
                                            method_cor='spearman',verbose=0){
    
    LS<-LayeredFromHierPath(hierarchy_name, datasets_path,dataset_name,
                            base_expr_data_path, base_expr_data, #one of those must not be null
                            base_has_decision, base_expr_RData,
                            method_cor,verbose)
   LS$RGM<- RandomGreyModules(LS$GrayArea, LS$expr_data, neg_cor_prop,
                             MAX_COR, MIN_COR, MaxSubmoduleSize, corpower,
                             sizeProbs) 
    prefix_path=paste0(datasets_path,'/',dataset_name,'/simulations/GraySubmodules/')
    dir.create(paste0(prefix_path, GraySubmodulesName), showWarnings=FALSE, recursive=TRUE)
    GraySpecs<- list(NAME=GraySubmodulesName,
                     base_dataset=dataset_name,
                     base_hierarchy=hierarchy_name,
                     neg_cor_prop=neg_cor_prop,
                     MAX_COR=MAX_COR,
                     MIN_COR=MIN_COR,
                     MaxSubmoduleSize=MaxSubmoduleSize,
                     corpower=corpower,
                     sizeProbs=sizeProbs)
    saveRDS(GraySpecs,paste0(prefix_path,GraySubmodulesName,'/Specs.rds'))
    return(LS)
}

ApplyGreyModulesSpecs<- function(layered_sim, GreySpecs)
{
    layered_sim$RGM<- RandomGreyModules(layered_sim$GrayArea, layered_sim$expr_data, GreySpecs$neg_cor_prop,
                                        GreySpecs$MAX_COR, GreySpecs$MIN_COR, GreySpecs$MaxSubmoduleSize,
                                        GreySpecs$corpower, GreySpecs$sizeProbs)
    return(layered_sim)
}

LoadGreyModulesSpecs<- function(GraySubmodulesName, dataset_name, datasets_path) {
    prefix_path=paste0(datasets_path,'/',dataset_name,'/simulations/GraySubmodules/')
    NewSpecs<-readRDS(paste0(prefix_path,GraySubmodulesName,'/Specs.rds'))
    return(NewSpecs)
}

ApplyGreyModulesSpecsFromFile<- function(layered_sim, GraySubmodulesName, dataset_name, datasets_path) {
   GreySpecs<- LoadGreyModulesSpecs(GraySubmodulesName, dataset_name, datasets_path)
   applied<- ApplyGreyModulesSpecs(layered_sim, GreySpecs)
   return(applied) 
}

EvaluateGreyModSpec<- function(GreySpec=NULL, GreySpecName=NULL, #one must not be null,
                               save_evaluation=FALSE, #if true GreySpecName must be not null
                               calcClCoef=FALSE, calcDTOM=FALSE,
                               layered_sim=NULL, hier=NULL, hierarchy_name=NULL, #one must not be null,
                               dataset_name=NULL, datasets_path=NULL, #if GreySpecFile or hierPath
                               base_expr_data_path=NULL, base_expr_data=NULL, #one of those must not be null if layered_sim=NULL
                               base_has_decision=FALSE, base_expr_RData=TRUE,
                               method_cor='spearman',verbose=0){
    if(is.null(GreySpec)) if(is.null(GreySpecName)) stop('one of GreySpec or GreySpecName must not be NULL')
                          else  GreySpec<- LoadGreyModulesSpecs(GreySpecName, dataset_name, datasets_path)
    if(is.null(layered_sim)) if(is.null(hier)){ if(is.null(hierarchy_name)) stop('One of: layered_sim, hier, hierarchy_name must not be null')
                                               else layered_sim<- LayeredFromHierPath(hierarchy_name, datasets_path,dataset_name,
                                                                   base_expr_data_path, base_expr_data, 
                                                                   base_has_decision, base_expr_RData,
                                                                   method_cor,verbose)
                                              }
                             else
                             layered_sim<- LayeredFromHierarchy(hier, datasets_path,
                                                                   base_expr_data_path, base_expr_data, 
                                                                   base_has_decision, base_expr_RData,
                                                                   method_cor,verbose)
    LS<-ApplyGreyModulesSpecs(layered_sim, GreySpec)
    if (is.null(hier)) hier<- ReadHierarchy(hierarchy_name, datasets_path, dataset_name)
    real_net_name<-hier$base$base_sim_specs$base_network_name
    tempc<-cor(LS$RGM, method=method_cor)
    tempadj<- abs(tempc)^strtoi(real_net_name)
    rm(tempc)
    GrayStats<- list(deg_distrGrey<- colSums(tempadj[,LS$GrayArea]))    
    if (calcClCoef){ print('Clustering coef...'); GrayStats$ClCoef<- clusterCoef(tempadj)[LS$GrayArea]}
    if (calcDTOM) { tempdTOM<-1-TOMsimilarity(tempadj);GrayStats$TOMdDistr<-  as.vector(tempdTOM[LS$GrayArea,LS$GrayArea] )  }
    if (save_evaluation){
        prefix_path=paste0(datasets_path,'/',dataset_name,'/simulations/GraySubmodules/')
        saveRDS(GrayStats,paste0(prefix_path,GreySpecName,'/Stats.rds'))}
    return(GrayStats)
}

LoadGreySpecStats<- function(GreySpecName, datasets_path, dataset_name)
{
    prefix_path=paste0(datasets_path,'/',dataset_name,'/simulations/GraySubmodules/')
        GrayStats = readRDS(paste0(prefix_path,GreySpecName,'/Stats.rds'))
        return(GrayStats)
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
    print("doing")
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


RandomSmallSubModulesFromHierPath<- function(hierarchy_name, datasets_path,dataset_name,
                                        MAX_COR,MIN_COR,MaxSubmoduleSize, neg_cor_prop,
                            save_specs=FALSE, NewSpecName=NULL, #if save name must not be null
                            corpower=1, sizeProbs=rep( 1/MaxSubmoduleSize,MaxSubmoduleSize),
                            GreySpecs=NULL, GreySpecsName=NULL,
                            base_expr_data_path=NULL, base_expr_data=NULL, #one of those must not be null
                           base_has_decision=FALSE, base_expr_RData=TRUE,
                           method_cor='spearman',verbose=0){

    LS<-LayeredFromHierPath(hierarchy_name, datasets_path,dataset_name,
                            base_expr_data_path, base_expr_data, #one of those must not be null
                            base_has_decision, base_expr_RData,
                            method_cor,verbose)
   print(!is.null( GreySpecs ))
   print(!is.null( GreySpecsName  ))
   print(!((is.null( GreySpecs )) & (is.null( GreySpecsName ))))
   if(!is.null( GreySpecs ) )LS<- ApplyGreyModulesSpecs (LS, GreySpecs)
   if(!is.null( GreySpecsName ) )LS<- ApplyGreyModulesSpecsFromFile (LS, GreySpecsName, dataset_name, datasets_path)
    if(!((is.null( GreySpecs )) & (is.null( GreySpecsName )))) {LS$SubMods<-  RandomSmallSubModules(LS$hierarchy$base$color_labels, LS$RGM, neg_cor_prop,
                                 MAX_COR, MIN_COR, MaxSubmoduleSize,
                                 sizeProbs=rep( 1/MaxSubmoduleSize,MaxSubmoduleSize)); print("not null opt")
   } else{print("null opt");
    LS$SubMods<-RandomSmallSubModules(LS$hierarchy$base$color_labels, LS$expr_data, neg_cor_prop,
                                 MAX_COR, MIN_COR, MaxSubmoduleSize,
                                 sizeProbs=rep( 1/MaxSubmoduleSize,MaxSubmoduleSize))}
    if (save_specs){
    RSpecs<- list(NAME=NewSpecName,
                     base_dataset=dataset_name,
                     base_hierarchy=hierarchy_name,
                     neg_cor_prop=neg_cor_prop,
                     MAX_COR=MAX_COR,
                     MIN_COR=MIN_COR,
                     MaxSubmoduleSize=MaxSubmoduleSize,
                     corpower=corpower,
                     sizeProbs=sizeProbs)
    prefix_path=paste0(datasets_path,'/',dataset_name,'/simulations/RandomProperSubmodules/')
    dir.create(paste0(prefix_path, NewSpecName), showWarnings=FALSE, recursive=TRUE)
    saveRDS(RSpecs,paste0(prefix_path,NewSpecName,'/Specs.rds'))
                    }
    return(LS)
}


ApplyRandModulesSpecs<- function(layered_sim, RSpecs, to_RGM=TRUE)
{
    if (to_RGM)layered_sim$SubMods<- RandomSmallSubModules(layered_sim$base$color_labels, layered_sim$RGM, RSpecs$neg_cor_prop,
                                 RSpecs$MAX_COR, RSpecs$MIN_COR, RSpecs$MaxSubmoduleSize,
                                 RSpecs$sizeProbs)
    else
        layered_sim$SubMods<- RandomSmallSubModules(layered_sim$base$color_labels, layered_sim$expr_data, RSpecs$neg_cor_prop,
                                 RSpecs$MAX_COR, RSpecs$MIN_COR, RSpecs$MaxSubmoduleSize,
                                 RSpecs$sizeProbs)
    return(layered_sim)
}

LoadRandomModulesSpecs<- function(RSubmodulesName, dataset_name, datasets_path) {
    prefix_path=paste0(datasets_path,'/',dataset_name,'/simulations/RandomProperSubmodules/')
    NewSpecs<-readRDS(paste0(prefix_path,RSubmodulesName,'/Specs.rds'))
    return(NewSpecs)
}

ApplyRandModulesSpecsFromFile<- function(layered_sim, RSubmodulesName, dataset_name, datasets_path,to_RGM=TRUE) {
   NewSpecs<- LoadRandomModulesSpecs(RSubmodulesName, dataset_name, datasets_path)
   applied<- ApplyRandModulesSpecs(layered_sim, NewSpecs,to_RGM)
   return(applied) 
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


CompareClusteringAndLayered<- function(clust_result,layered_sim,
side="simulation", main="sim vs base clusterings",...)
{
plotDendroAndColors(clust_result$geneTree, 
as.data.frame(c(list(clust_result$color_labels),
  lapply(layered_sim$hierarchy, function(x) x$color_labels))),c(side,names(layered_sim$hierarchy)), dendroLabels=FALSE, main=main, ...)

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
