

dissTOM2dendro <- function (distance_matrix, method="average"){
hclust(distance_matrix, method = method)
}

dendro2clusters <- function (geneTree, distance_matrix, sensivity = 2, pam2parent= FALSE,
                             minModuleSize=30, doPAM=TRUE, verbosity_level=0, ...) {
  cutreeDynamic(dendro = geneTree, distM = distance_matrix,
                deepSplit = sensivity, pamRespectsDendro = pam2parent,
                minClusterSize = minModuleSize, 
                pamStage= doPAM, verbose = verbosity_level, ...)
}

colors2MEs <- function (clabels, dE, ...) {
  moduleEigengenes(dE, clabels, ...)
}

#Starts with adjacency matrix (i.e. output from adjacency.fromSimilarity from WGCNA)
#calculation of TOM in case of diss_as_TOM=TRUE is compute-intensive
#if that's the problem consider using ClusteringResults.fromSimilarity
ClusteringResults.fromAdjacency <- function ( adj_mat, dataExpr, save_steps=FALSE, save_final=TRUE,
                                             diss_as_TOM= TRUE,ts_other_args=NULL,
                                             hclust_method = "average", minModSize=30,
                                              dc_sensivity= 2, dc_pam2parent=FALSE, dc_doPAM = TRUE, 
                                              verbose_lvl=2, dc_other_args=NULL, ME_other_args= NULL,
                                             calculateMEs=TRUE){
  
  
  if (verbose_lvl>0) print("Calculating dissimilarity...")
  if (diss_as_TOM){
  TOM_args <- list( adj_mat)
  if (!is.null(ts_other_args)) TOM_args <- c(TOM_args, ts_other_args)
  TOM_args <- c(TOM_args, list(verbose= verbose_lvl))
  if (verbose_lvl==4) {print("Supplied TOMsimilarity args: "); print(TOM_args[-1]) }
  
  TOM <- do.call(TOMsimilarity, TOM_args)
  } else TOM <- adj_mat
  dissTOM <- 1-TOM
  
  if (verbose_lvl>0) print("Calculated dissimilarity")
  
  if (save_steps) {if (verbose_lvl>0) print("Saving dissimilarity matrix..."); saveRDS(dissTOM, "dissMat.rds");
                   if (verbose_lvl>0) print("Saved dissimilarity matrix to dissMat.rds")}
  
  if (verbose_lvl>0) print("Making dendogram...")
  dendogram <- dissTOM2dendro(as.dist(dissTOM), method= hclust_method)
  
  if (verbose_lvl>0) print("Finished dendogram.")
  
  if (save_steps) {if (verbose_lvl>0) print("Saving dendrogram..."); saveRDS(dendogram, "geneTree.rds");
    if (verbose_lvl>0) print("Saved dendogram to geneTree.rds")}
  
  if (verbose_lvl>0) print("Calculating cluster labels...")
  d2c_args <- list(dendogram, dissTOM, 
                   sensivity=dc_sensivity, pam2parent= dc_pam2parent,
                   minModuleSize=minModSize, doPAM=dc_doPAM, verbosity_level=verbose_lvl )
  if(!is.null(dc_other_args)) d2c_args <- c(d2c_args, dc_other_args)
  if (verbose_lvl==4) {print("Supplied cutreeDynamic args: "); print( d2c_args[-c(1,2)]) }
  
  dynamicMods <- do.call(dendro2clusters, d2c_args)
  
  if (verbose_lvl>0) print("Finished clustering.")
  
  if (save_steps) {if (verbose_lvl>0) print("Saving cluster labels..."); saveRDS(dynamicMods, "dynamicMods.rds");
    if (verbose_lvl>0) print("Saved cluster labels to dynamicMods.rds")}
  
  
  dynamicColors <- labels2colors(dynamicMods)
  
  if (calculateMEs){
  if (verbose_lvl>0) print("Calculating module eigengenes.")
  c2M_args= list( dynamicColors, as.data.frame(dataExpr), verbose= verbose_lvl)
  if(!is.null(ME_other_args)) c2M_args <- c(c2M_args, ME_other_args)
  if (verbose_lvl==4) {print("Supplied moduleEigengenes args: "); print( c2M_args[-c(1,2)]) }
  
  ME_list <- do.call(colors2MEs, c2M_args )
  
  if (verbose_lvl>0) print("Calculated Eigengenes.")
  
  if (save_steps) {if (verbose_lvl>0) print("Saving module eigengene data..."); saveRDS(dynamicMods, "ME_list.rds");
    if (verbose_lvl>0) print("Saved module eigengene data to ME_list.rds")}
  } else{
  ME_list <- "X"}
  
  final_results <- list(geneTree = dendogram,
                        color_labels = dynamicColors,
                        ME_data = ME_list
                        )
  
  if (save_final) {if (verbose_lvl>0) print("Saving final results (dendogram, color labels, MEs..."); saveRDS(final_results, "Adj2Clust_results.rds");
    if (verbose_lvl>0) print("Saved final results to Adj2Clust_results.rds")}
  
  return(final_results)
}

#Here similarity is meant as a measure of clustering

ClusteringResults.fromSimilarity<- function(sim_mat, dataExpr, save_steps=FALSE, save_final=TRUE,
                                 hclust_method = "average", minModSize=30,
                                 dc_sensivity= 2, dc_pam2parent=FALSE, dc_doPAM = TRUE, 
                                 verbose_lvl=2, dc_other_args=NULL, ME_other_args= NULL){
  dissTOM <- 1-sim_mat
  
  
  if (verbose_lvl>0) print("Making dendogram...")
  dendogram <- dissTOM2dendro(as.dist(dissTOM), method= hclust_method)
  
  if (verbose_lvl>0) print("Finished dendogram.")
  
  if (save_steps) {if (verbose_lvl>0) print("Saving dendrogram..."); saveRDS(dendogram, "geneTree.rds");
    if (verbose_lvl>0) print("Saved dendogram to geneTree.rds")}
  
  if (verbose_lvl>0) print("Calculating cluster labels...")
  d2c_args <- list(dendogram, dissTOM, 
                   sensivity=dc_sensivity, pam2parent= dc_pam2parent,
                   minModuleSize=minModSize, doPAM=dc_doPAM, verbosity_level=verbose_lvl )
  if(!is.null(dc_other_args)) d2c_args <- c(d2c_args, dc_other_args)
  if (verbose_lvl==4) {print("Supplied cutreeDynamic args: "); print( d2c_args[-c(1,2)]) }
  
  dynamicMods <- do.call(dendro2clusters, d2c_args)
  
  if (verbose_lvl>0) print("Finished clustering.")
  
  if (save_steps) {if (verbose_lvl>0) print("Saving cluster labels..."); saveRDS(dynamicMods, "dynamicMods.rds");
    if (verbose_lvl>0) print("Saved cluster labels to dynamicMods.rds")}
  
  
  dynamicColors <- labels2colors(dynamicMods)
  
  
  if (verbose_lvl>0) print("Calculating module eigengenes.")
  c2M_args= list( dynamicColors, as.data.frame(dataExpr), verbose= verbose_lvl)
  if(!is.null(ME_other_args)) c2M_args <- c(c2M_args, ME_other_args)
  if (verbose_lvl==4) {print("Supplied moduleEigengenes args: "); print( c2M_args[-c(1,2)]) }
  
  ME_list <- do.call(colors2MEs, c2M_args )
  
  if (verbose_lvl>0) print("Calculated Eigengenes.")
  
  if (save_steps) {if (verbose_lvl>0) print("Saving module eigengene data..."); saveRDS(dynamicMods, "ME_list.rds");
    if (verbose_lvl>0) print("Saved module eigengene data to ME_list.rds")}
  
  
  final_results <- list(geneTree = dendogram,
                        color_labels = dynamicColors,
                        ME_data = ME_list
  )
  
  if (save_final) {if (verbose_lvl>0) print("Saving final results (dendogram, color labels, MEs..."); saveRDS(final_results, "Adj2Clust_results.rds");
    if (verbose_lvl>0) print("Saved final results to Adj2Clust_results.rds")}
  
  return(final_results)
}
  
