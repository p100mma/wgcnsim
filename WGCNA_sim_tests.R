


Repeat_simulation_from_input <- function(times=12, input_args_list, verbose=0, other_named_args=NULL) {
  sim_list <- vector( mode = "list", length= times)
  for (i in c(1:times)){
    sim_list[[i]] <- simulateDatExpr_fromInput(input_args_list,verbose, other_named_args)
  }
  return(sim_list)
}

reorder_sim_results <- function(sim_results, real_clustering_result) {
  lapply(sim_results, ReorderSimByReal, real_clustering_result)
}

#RAM friendly functions 

compute_and_save_cors <- function(ordered_simExprs, method_cor, verbose=0, name_prefix="sim_", ...){
  i=1
  save_names <- paste0(name_prefix,c(1:length(ordered_simExprs)),"_",method_cor,".rds")
  for (ose in ordered_simExprs){
    if (verbose>0) print(i)
    cur_cor<- cor(ose, method=method_cor, ...)
    saveRDS(cur_cor,save_names[[i]])
    rm(cur_cor)
    i=i+1
  }
  invisible( save_names)
}

load_cors_compute_save_adj<- function(load_names, adj_power=5, verbose=0, name_prefix="sim_", ...){
  save_names<- paste0(name_prefix,c(1:length(load_names)),"_adjacency.rds")
  for (i in seq_along(load_names)){
    if (verbose>0) print(i)
    cur_cor <- readRDS(load_names[[i]])
    cur_adj <- adjacency.fromSimilarity(cur_cor, power=adj_power, ...)
    saveRDS(cur_adj,save_names[[i]])
    rm(cur_cor)
    rm(cur_adj)
  }
  invisible(save_names)
}

load_adjs_compute_save_TOMs<- function(load_names, verbose=0,name_prefix="sim_", ...) {
  save_names<- paste0(name_prefix,c(1:length(load_names)),"_TOM.rds")
  for (i in seq_along(load_names)){
    if (verbose>0) print(i)
    cur_adj <- readRDS(load_names[[i]])
    cur_TOM <- TOMsimilarity(cur_adj, ...)
    saveRDS(cur_TOM,save_names[[i]])
    rm(cur_adj)
    rm(cur_TOM)
  }
  invisible(save_names)
}

load_TOMs_compute_colors <- function(load_names, hclust_method = "average", minModSize=86,
                                     dc_sensivity= 2, dc_pam2parent=FALSE, dc_doPAM = TRUE, 
                                     verbose_lvl=2, dc_other_args=NULL) {
  sim_colors <- vector( mode="list", length=length(load_names))
  for (i in seq_along(load_names)){
    if (verbose_lvl>0) print(i)
    cur_TOM <- readRDS(load_names[[i]])
    cur_dissTOM <- 1- cur_TOM
    cur_dendro <- dissTOM2dendro(as.dist(cur_dissTOM), method= hclust_method)

    if (verbose_lvl>0) print("Calculating cluster labels...")
    d2c_args <- list(cur_dendro, cur_dissTOM, 
                     sensivity=dc_sensivity, pam2parent= dc_pam2parent,
                     minModuleSize=minModSize, doPAM=dc_doPAM, verbosity_level=verbose_lvl )
    if(!is.null(dc_other_args)) d2c_args <- c(d2c_args, dc_other_args)
    if (verbose_lvl==4) {print("Supplied cutreeDynamic args: "); print( d2c_args[-c(1,2)]) }
    cur_mods <- do.call(dendro2clusters, d2c_args)
    sim_colors[[i]]<- labels2colors(cur_mods)
    rm(cur_TOM)
    rm(cur_dissTOM)
    rm(cur_dendro)
    rm(cur_mods)
  }
  return(sim_colors)
}

Pairwise_labels_comparsion <- function(label_vectors, what_fun){
  idx_pair_apply <- function (x, y, what_fun, what) mapply( what_fun, what[x], what[y])
  idxs<- seq_along(label_vectors)
  matrix_result <- outer(idxs,idxs, idx_pair_apply, what_fun, label_vectors)
  rownames(matrix_result)<- names(label_vectors)
  colnames(matrix_result)<- names(label_vectors)
  return(matrix_result)
}


Resample_rows <- function(datExpr, n=NULL) {
  if (is.null(n)) n<- nrow(datExpr)
  sampled_idx <- sample(c(1:n),n, replace=TRUE)
  datExpr[sampled_idx,]
}

#do all of above but dont save intermediate steps on hard drive. Do each cor > adj > TOM sequentially
#saves final sim_expr_data (reordered) on hard drive and their clustering results in one file, one file for each repeat
TestRepeatedSim.RAMfriendly<- function(input_args_list, other_named_args, real_clustering_result,
                                       method_cor="spearman", adj_power=5,
                                       diss_as_TOM= TRUE,ts_other_args=NULL,
                                       hclust_method = "average", minModSize=30,
                                       dc_sensivity= 2, dc_pam2parent=FALSE, dc_doPAM = TRUE, 
                                       dc_other_args=NULL, ME_other_args= NULL,
                                       verbose=0,
                                       times=12,
                                       name_prefix="sim_and_clres") {  
  for (j in c(1:times)){
    if (verbose >0) print(paste0("Starting simulation number ",j))
    jth_cluster<- fromInput_simulate2ClusteringResults(input_args_list,other_named_args, real_clustering_result,
                                                    method_cor, adj_power, 
                                                    FALSE, FALSE,
                                                    diss_as_TOM,ts_other_args,
                                                    hclust_method, minModSize,
                                                    dc_sensivity, dc_pam2parent, dc_doPAM, 
                                                    dc_other_args, ME_other_args,
                                                    verbose=0)
    saveRDS(jth_cluster, paste0(name_prefix,j))
    rm(jth_cluster)
  }
}

#Will most likely overload memory on local PC

TestRepeatedSim<- function(real_clustering_result, simulation_list, cor_method="spearman", adj_power=5, clst_adj_other_args=NULL,
                            verbose=0){
  ordered_simulations_list<- lapply(simulation_list, ReorderSimByReal, real_clustering_result)
  sim_dataExpr_list<- lapply(ordered_simulations_list, function (x) x$datExpr )
  #sim_dataCors_list<- lapply(sim_dataExpr_list, cor, method=cor_method)
  #sim_dataAdjs_list<- lapply(sim_dataCors_list, adjacency.fromSimilarity, power= adj_power)
  sim_dataClst_list<- vector( mode="list", length = length(sim_dataExpr_list))
  
  if (verbose>0) print("Working on clustering simulations...")
  for (k in seq_along(sim_dataExpr_list)){
    if (verbose>1) print(paste0("Starting cluster ",k))
    print(paste("ncol:",ncol(cor(sim_dataExpr_list[[k]]))))
    cur_cor<- cor(sim_dataExpr_list[[k]], method=cor_method)
    gc()
    if (verbose>2) print(paste0("Done correlation ",k))
    cur_adj<- adjacency.fromSimilarity(cur_cor, power=adj_power)
    gc()
    if (verbose>3) print(paste0("Done adjacency ",k))
    cl_args <- list( cur_adj, sim_dataExpr_list[[k]])
    rm(cur_cor)
    rm(cur_adj)
    cl_args<- c(cl_args, list(verbose_lvl=verbose), clst_adj_other_args)
    sim_dataClst_list[[k]] <- do.call(ClusteringResults.fromAdjacency, cl_args)
    gc()
  }
  if (verbose>0) print("Done clustering simulations.")
  
  sim_dataLabs_list<- lapply(sim_dataClst_list, function(x) x$color_labels)
  Real_Sim_labels_list <- c( list(real_clustering_result$color_labels), sim_dataLabs_list )
  names(Real_Sim_labels_list) <- c("Real", paste0("simulation_",c(1:length(sim_dataLabs_list))))
  idx_pair_apply <- function (x, y, what_fun, what) mapply( what_fun, what[x], what[y])
  Labs_idx<- seq_along(Real_Sim_labels_list)
  RealSimRandMat<- outer(Labs_idx,Labs_idx, idx_pair_apply, RealSimRandIndex, Real_Sim_labels_list)
  rownames(RealSimRandMat)<- names(Real_Sim_labels_list)
  colnames(RealSimRandMat)<- names(Real_Sim_labels_list)
  repeat_results<- list(sim_cluster_results= sim_dataClst_list, 
                        sim_Exprdata_list= sim_dataExpr_list,
                        sim_dataLabs_list= sim_dataLabs_list,
                        RealSimRandIndexMat = RealSimRandMat)
  return(repeat_results)
}

##Functions for parameter searching

# #apply agg_func on outputs of fundamentalNetworkConpcepts that have length >1
# #returns numeric named vector
# fundamental.aggregate<- function( fundamental_stats, agg_func=mean) {
#   c(agg.Connectivity= agg_func(fundamental_stats$Connectivity),
#     agg.ScaledConnectivity = agg_func(fundamental_stats$ScaledConnectivity),
#     agg.ClusterCoef = agg_func(fundamental_stats$ClusterCoef),
#     agg.MAR = agg_func(fundamental_stats$MAR),
#     Density=fundamental_stats$Density,
#     Centralization=fundamental_stats$Centralization,
#     Heterogeneity=fundamental_stats$Heterogeneity)
# }
# 
# 
# #input: output from FundamentalNetworkConceps in form of numeric vector,
# #           where elements of original output that have length >1 are replaced by aggregate(i.e. mean),
# #       adjacency matrix of dataset to compare with
# #output: difference in aggregate parameters (from FundamentalNetworkConcepts in WGCNA)
# compute_param_difference<- function( reference_agg.stats, cur_adj, agg_func=mean){
#   cur.fundamental <- fundamentalNetworkConcepts(cur_adj)
#   cur.fundamental.agg <- fundamental.aggregate(cur.fundamental, agg_func)
#   reference_agg.stats- cur.fundamental.agg
# }
# 
# 
# 
# #returns aggregated fundamentalNetworkConcepts of simulation produced with supplied args,
# #and its Rand Index in relation with real_clustering_result
# #ref. prefixed arguments relate to original data
# #method_cor_construct determines the cor method argument used in estimating min and max cor, args to 
# #    simulateDatExpr
# #method_cor_clust determines the cor method argument used for construction of adjacency... etc. to cluster
# #    simulated data
# #customfunc is a logical argument, if TRUE, simulate data by using simulateDatExpr.customfunc,
# #if customfunc is TRUE, arguments eigen_relation and relation_params should be specified for it to work
# #as intended
# # ref.agg.stats is result of fundamental.aggregate applied to real data
# #agg_func is the function used in fundamental.aggregate on simulated data
# #adj_power is power used in constructing adjacency matrix from correlation matrix of simulated data
# #other arguments are named just like the origina arguments in functions:
# #  clust_res2simDEargs, simulateDatExpr_fromInput, ClusteringResults.fromAdjacency 
# 
# param_search.inputs2metrics<- function( ref.datExpr, ref.color_labels, ref.cormat, ref.MEs,
#                                         min_from_perm, abs_cor,
#                                         true_grey_frac, method_cor_construct,
#                                         other_named_args, customfunc, 
#                                         real_clustering_result, method_cor_clust,
#                                         adj_power, ref.agg.stats, agg_func,
#                                         diss_as_TOM,ts_other_args,
#                                         hclust_method, minModSize,
#                                         dc_sensivity, dc_pam2parent, dc_doPAM, 
#                                         dc_other_args, ME_other_args,
#                                         eigen_relation = '+', relation_params=NULL){
#   sim_args<- clust_res2simDEargs(ref.datExpr, ref.color_labels, ref.cormat, ref.MEs, 
#                                  min_from_perm, abs_cor,
#                                  true_grey_frac, 
#                                  FALSE,TRUE, 0,  method=method_cor_construct)
#   if (!customfunc)
#   sim_res<-simulateDatExpr_fromInput(sim_args, 0, other_named_args)
#   else
#   sim_res<-simulateDatExpr_fromInput.customfunc(sim_args, eigen_relation, relation_params, 0, other_named_args)
#   rm(sim_args)
#   sim_res <- sim_res$sim_result
#   order_sim_res <- ReorderSimByReal(sim_res, real_clustering_result)
#   sim_res_expr_data <- order_sim_res$datExpr
#   rm(order_sim_res)
#   n_nan<- sum(is.na(sim_res_expr_data))
#   if (n_nan >0 ){ which_nan <- is.na(sim_res_expr_data); nanlessmax<-max( sim_res_expr_data[!which_nan]);
#                   sim_res_expr_data[which_nan] <- nanlessmax}
#   sim_cor <- cor(sim_res_expr_data, method=method_cor_construct )
#   sim_adj <- adjacency.fromSimilarity(sim_cor, power=adj_power); rm(sim_cor)
#   fundamental.score<- compute_param_difference(ref.agg.stats, sim_adj, agg_func)
#   sim_cluster <- ClusteringResults.fromAdjacency( sim_adj, ref.datExpr, FALSE, FALSE,
#                                                                diss_as_TOM,ts_other_args,
#                                                                hclust_method, minModSize,
#                                                                dc_sensivity, dc_pam2parent, dc_doPAM, 
#                                                                verbose_lvl=0, dc_other_args, ME_other_args)
#   RI <- RealSimRandIndex(real_clustering_result$color_labels, sim_cluster$color_labels  );
#   rm(sim_cluster); rm(sim_res_expr_data)
#   c(as.list(fundamental.score),
#        list(RI=RI, nan_replacements=n_nan))
# }
# #try every combination of parameters in params_list and compute fundamentalNetworkConcepts
# #associated with said combination and RI, in relation to real data
# #ref. prefixed arguments relate to original data
# #real_clustering_result is result of ClusteringResults functions 
# #ref.agg.stats is result of fundamental.aggregate applied to real data
# #params_list is a list of lists,
# #should contain lists of named elements with one element 
# #for every argument to param_search.inputs2metrics,
# #     apart from those already present in function definition
# #     every sublist should be the same length
# #     list of lists is basically cartesian product of vectors of parameter values
# #     expand.grid can be used to construct the lists.
# 
# param_search.simulation<- function(ref.datExpr, ref.color_labels, ref.cormat, ref.MEs,
#                         real_clustering_result, ref.agg.stats, 
#                         params_list ){
#     results_list <- vector( mode= "list", length= length(params_list))
#     i=1
#   for (params in params_list){
#     print(paste0("Working on ",i))
# cur_metrics <- param_search.inputs2metrics( ref.datExpr, ref.color_labels, ref.cormat, ref.MEs,
#                                             params$min_from_perm, params$abs_cor,
#                                             params$true_grey_frac, params$method_cor_construct,
#                                             params$other_named_args, params$customfunc, 
#                                             real_clustering_result, params$method_cor_clust,
#                                             params$adj_power,ref.agg.stats, params$agg_func,
#                                             params$diss_as_TOM,params$ts_other_args,
#                                             params$hclust_method, params$minModSize,
#                                             params$dc_sensivity, params$dc_pam2parent, params$dc_doPAM, 
#                                             params$dc_other_args, params$ME_other_args,
#                                             params$eigen_relation, params$relation_params)
# ith_result <- c(params, cur_metrics)
# results_list[[i]]<- ith_result
#   i=i+1
#   }
#     return(results_list)
# }
# 






