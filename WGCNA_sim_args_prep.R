

 ncol_check<- function(left, right, left_name=NULL, right_name=NULL){
   message=NULL
   OK <- (ncol(left)== ncol(right))
   if (is.null(left_name) || is.null(right_name) )
   {
   if (!OK) message <- paste0("Error, left has ",
                              ncol(left)," columns while right has ",
                              ncol(right)," columns.")}
   else{
     if (!OK)  message <- paste0("Error, ",left_name," has ",
                                ncol(left)," columns while ",right_name," has ",
                                ncol(right)," columns.")}
   list( status = OK, msg=message )
 }
 
 ncol_len_check<- function(left, right, left_name=NULL, right_name=NULL){
   message=NULL
   OK <- (ncol(left)== length(right))
   if (is.null(left_name) || is.null(right_name) )
   {
     if (!OK) message <- paste0("Error, left has ",
                                ncol(left)," columns while right has length of ",
                                length(right))}
   else{
     if (!OK)  message <- paste0("Error, ",left_name," has ",
                                 ncol(left)," columns while ",right_name," has length of ",
                                 length(right))}
   list( status = OK, msg=message )
 }

cor2sim_neg_cor_prop <- function(cormat){
  sum(cormat<0)/(ncol(cormat)*ncol(cormat))
}

ME_module_gene_cors<- function(sDE, ME, method="pearson"){
  #calculates a vector of correlations between every column of sDE (multiple columns) and ME (1 column vector)
  unlist(lapply(sDE, cor, ME, method=method))
}

matfunc_on_module<- function(colorname, dEdf, clabels, matfunc, MEs=NULL, verbose=0,...){
  #based on labeling made by clustering and labels2colors,
  #select genes labelled with 'colorname' from 'dEdf' (Expression data)
  #and do 'matfunc' on them. 'MEs' are module eigengenes- matfunc can also work with
  #supplied eigengene, identified by 'colorname'
  #dEdf should be list of vectors (or data frame)
  
  err_check<- ncol_len_check(dEdf, clabels, "expression data", "module labels vector")
  if (!err_check$status) stop(err_check$msg)
  
  if (!is.null(MEs)) cMEname <- paste0("ME",colorname)
  idx <- clabels==colorname
  subdE <- dEdf[,idx]
  if (is.null(MEs)) res<- matfunc(subdE, ...) else res<- matfunc(subdE, MEs[[cMEname]], ...)
  if (verbose>0) print(paste0("Done calculating for ",colorname))
  return(res)
}

colors2ME_gene_cors <- function(colors2apply, dE, clabels, Eigengenes, verbose=0, ...){
  #get correlations of Eigengene with its module genes for each color of colors2apply
  #dE is expression data, clabels are labels of modules in whole dataset (length(clabels)==ncol(dE))
  #further arguments (...) are delivered to cor (i.e. method="spearman")
  #returns list of vectors, each vector has length equal to size of module from colors2apply
 lapply(colors2apply, matfunc_on_module, as.data.frame( dE), clabels, ME_module_gene_cors,Eigengenes, verbose, ...)
}

permute_column_wise <- function(dE){
  #permute values in columns of dE independently and return altered data
  for (i in c(1:ncol(dE))){
    idx <- sample(c(1:nrow(dE)), nrow(dE))
    dE[,i] = dE[idx,i]
  }
  return(dE)
}

colors2sim_cor_params <- function(colors2apply, datExpr, clabels, Eigengenes, abs_cor=TRUE, min_from_perm=FALSE, verbose=0, ...){
  #get min and max of correlations from of Eigengene with its module genes for each color of colors2apply
  #datExpr is expression data, clabels are labels of modules in whole dataset (length(clabels)==ncol(datExpr))
  #further arguments (...) are delivered to cor (i.e. method="spearman")
  #returns list of vectors of same length which is equal to length of colors2apply
  #if min_from_perm is TRUE, minimal correlation in each module of colors2apply is the largest of 
  #     {90th quantile of cors with Eigengene in module with column values permuted independently,
  #      real estimated minimal correlation obtained from calculations}
  #if abs_cor==TRUE, then returned min and max relates to absolute values of correlations
  corvecs <- colors2ME_gene_cors(colors2apply, datExpr, clabels, Eigengenes,verbose, ...)
  if (abs_cor) corvecs <- lapply(corvecs,abs)
  min_<- unlist(lapply(corvecs, min))
  max_<- unlist(lapply(corvecs, max))
  if (min_from_perm){
    if (verbose>0) print("Calculating module-minimal correlations from column-wise permuted data")
    pdatExpr<-permute_column_wise(datExpr)
    rand_corvecs<-colors2ME_gene_cors(colors2apply, pdatExpr, clabels, Eigengenes, verbose, ...)
    if (abs_cor) rand_corvecs <- lapply(rand_corvecs, abs)
    max_rand_cor<-unlist(lapply(rand_corvecs, quantile, 0.9))
    min_<-ifelse(min_<max_rand_cor, max_rand_cor, min_ )
  }
  return(list(min_cor=min_, max_cor=max_))
}

colors2props <- function( clabels, distinct_ordered_colors=NULL,true_grey_frac=0.5){
 #get proportions of modules in a form suitable for simulateDatExpr,
 #from module labels produced by WGCNA package. clabels should contain label for every gene 
 #in original data.
 #true_grey_frac is proportion of "true grey genes" (read below)
 #simulateDatExpr requires that the last entry of proportions vector should be proprtion of 
 #"true grey genes", that is genes truly independent of any module. proportions should also
 #sum to quantity smaller than 1 (1 - sum(proportions) is proportion of genes near modules,
 #   weakly correlated to Eigengenes)
  module_props<- table(clabels)
  gray_prop<- module_props[["grey"]]
  module_props<- module_props[ names(module_props)!="grey"]
  module_props<- c( module_props, gray_prop)
  module_props[[length(module_props)]]<- module_props[[length(module_props)]]*true_grey_frac
  module_props<- module_props/length(clabels)
  if (!is.null(distinct_ordered_colors)){
    wout_grey<-module_props[c(1:(length(module_props)-1))]
    wout_grey<-wout_grey[order(names(wout_grey), distinct_ordered_colors )]
    module_props[c(1:(length(module_props)-1))]<-wout_grey
  }
  return(module_props)
}

clust_res2simDEargs <- function(dataExpr, color_labels, cormat, MEs=NULL, 
                                min_from_perm=FALSE, abs_cor=TRUE,
                                true_grey_frac=0.5, 
                                save_steps=FALSE,save_final=TRUE, verbose=0,
                                ME_data_path=NULL,
                                ME_dataInFolderTree=FALSE, datasets_path=NULL,
                                dataset_name=NULL, network_name=NULL, clustering_name=NULL ){
  #get the input arguments to the simulateDatExpr from 
  #expression data (dataExpr) and results of clustering by WGCNA package
  #MEs are eigengenes, if NULL, then either:
  #     - ME_data_path must have .rds file with ME_data (result of moduleEigengenes from WGCNA),
  #     - ME_dataInFolderTree must be TRUE, and arguments datasets_path,dataset_name,network_name,
  #       ... and clustering_name must lead to ME_data.rds, : 
  #           datasets_path/dataset_name/networks/network_name/clusterings/clustering_name/ME_data.rds
  #cormat is original correlation matrix from which clustering was derived
  #color_labels are labels of genes produced by labels2colors from WGCNA package
  #if min_from_perm is TRUE, estimated minimal correlation in each module is the largest of 
  #     {90th quantile of cors with Eigengene in module with column values permuted independently,
  #      real minimal correlation obtained from calculations}
  #true_grey_frac determines the proportion of "true grey genes", last entry of modProportions argument to...
  #   ...simulateDatExpr
  #if abs_cor==TRUE, then returned min and max relates to absolute values of correlations
  #
  #returns a list of arguments named as correspoding formals in simulateDatExpr
  #list contains vectors of minimal and maximal correlations, gene means, eigengenes,
  #proportions of modules (in terms of gene counts, without "grey" module as in simulateDatExpr modProportions),
  #and nGenes which is ncol(dataExpr). nGenes in simulateDatExpr could be different than calculated here,
  # but then length(geneMeans) must equal to nGenes, 
  # so in that case geneMeans must be modified accordingly in future call of simulateDatExpr
  err_check<- ncol_check(dataExpr, cormat, "expression data", "correlation matrix")
  if (!err_check$status) stop(err_check$msg)
  if (is.null(MEs)) if(!is.null(ME_data_path)) MEs<- readRDS(ME_data_path)$eigengenes
                    else { if(!ME_dataInFolderTree) stop(" if MEs=NULL, then either:
       - ME_data_path must have .rds file with ME_data (result of moduleEigengenes from WGCNA),
       - ME_dataInFolderTree must be TRUE, and arguments datasets_path,dataset_name,network_name,
         ... and clustering_name must lead to ME_data.rds, : 
             datasets_path/dataset_name/networks/network_name/clusterings/clustering_name/ME_data.rds")
                          prefix_path=paste0(datasets_path,'/',dataset_name,'/networks/',network_name,'/clusterings/',clustering_name,'/')
                          MEs<- readRDS(paste0(prefix_path,'ME_data.rds'))$eigengenes
                         }
  MEnames <- names(MEs)
  colornames<- sub(".*ME", "", MEnames)
  colornames<- colornames[ colornames!="grey"]     #skip grey module (independent genes)
  
  
  if (verbose>0) print("Calculating min and max correlations.")
  cor_par_list <- colors2sim_cor_params(colornames, dataExpr, color_labels, MEs,abs_cor,min_from_perm,verbose  )
  if (verbose>0) print("Done calculations of min and max correlations.")
  
  if(save_steps) {if (verbose>0) print("Saving min and max correlations..."); saveRDS(cor_par_list, "cor_pars.rds")
                  if (verbose>0) print("Saved to file cor_pars.rds")}
  
  
  if (verbose>0) print("Calculating the proportion of negative correlations.")
  neg_cor_prop <- cor2sim_neg_cor_prop(cormat)
  if (verbose>0) print("Done calculation of negative correlations prop.")
  
  if(save_steps) {if (verbose>0) print("Saving proportion of negative orrelations...")
                  saveRDS(neg_cor_prop, "neg_cor_prop.rds")
                  if (verbose>0) print("Saved to file neg_cor_prop.rds")}
  
  
  if (verbose>0) print("Calculating proportions of genes in modules.")
  mod_props<- colors2props(color_labels, colornames, true_grey_frac)
  if (verbose>0) print("Done calculation of proportions of genes in modules.")
  
  if(save_steps) {if (verbose>0) print("Saving proportions of modules...")
    saveRDS(mod_props, "mod_props.rds")
    if (verbose>0) print("Saved to file mod_props.rds")}
  
  
  if (verbose>0) print("Calculating gene means.")
  gmeans <- unlist(lapply(as.data.frame(dataExpr),mean))
  if (verbose>0) print("Done calculation of gene means.")
  
  if(save_steps) {if (verbose>0) print("Saving gene means...")
    saveRDS(gmeans, "gMeans.rds")
    if (verbose>0) print("Saved to file gMeans.rds")}
  
  
  final_args<- list(eigengenes= MEs[names(MEs)!="MEgrey"], nGenes=ncol(dataExpr), 
                    modProportions=mod_props,
                    minCor=cor_par_list$min_cor, maxCor=cor_par_list$max_cor,
                    propNegativeCor=neg_cor_prop,
                    geneMeans=gmeans)
  
  if (save_final) {if (verbose>0) print("Saving final results"); saveRDS(final_args, "simArgs.rds");
    if (verbose>0) print("Saved final results to simArgs.rds")}
  
  return(final_args)
  
}


  
