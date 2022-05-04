


len_check<- function(left, right, left_name=NULL, right_name=NULL){
  message=NULL
  OK <- (length(left)== length(right))
  if (is.null(left_name) || is.null(right_name) )
  {
    if (!OK) message <- paste0("Error, left has length of ",
                               length(left)," while right has length of ",
                               length(right))}
  else{
    if (!OK)  message <- paste0("Error, ",left_name," has length of ",
                                length(left)," while ",right_name," has length of ",
                                length(right))}
  list( status = OK, msg=message )
}


simulateDatExpr_fromInput<- function(input_args_list, verbose=0, other_named_args=NULL){
  #calls function simulating expression data from input_args_list computed beforehand
  #by clust_res2simDEargs (or prepared independently).
  #input_args_list contains arguments: eigengenes, nGenes, modProportions, minCor, maxCor, propNegativeCor, geneMeans
  #or less (all of them if directly using result of clust_res2simDEargs). remaining arguments excluding 'verbose'
  #are to be put in other_named_args
  #returns output of simulateDatExpr and list of supplied arguments (containing only those given and 'verbose')
  if (sum(names(input_args_list) %in% names(other_named_args))!=0) stop(
    "Error: common members in input_args_list and other_named_args")
  if (sum("verbose" %in% union(names(other_named_args), names(input_args_list)))!=0) stop(
    "Error: put 'verbose' seperately as argument, not as a member of input_args_list or other_named_args")
  args_<- c(input_args_list, list(verbose= verbose))
  args_<- c(args_, other_named_args)
  if (verbose>0) {print("Calling simulateDatExpr with given args: "); print(args_)}
  sim_result<-do.call(simulateDatExpr, args_)
  return(list(sim_result=sim_result, simulation_call_args=args_))
}

same_members_check <- function(A,B, A_name=NULL, B_name=NULL){
  #check if two vectors have same elements
  OK<-!( length(intersect(A,B))!=length(A) || length(intersect(A,B))!=length(B))
  msg <- NULL
  if (!OK){
    if (is.null(A_name)) msg<- paste0("Error: left has ", length(A)," members while right has ", length(B))
    else msg<-paste0("Error: ",A_name," has ", length(A)," members while ",B_name," has ",length(B))
  }
    list(status=OK, msg=msg )
    }

members_count_check <- function(sim_labels, real_labels, sim_name, real_name){
  #check if two vectors have same counts of each element
  #to be called after same_members_check
  real_labels_counts<-table(real_labels)
  sim_labels_counts<-table(sim_labels)
  err_len<- len_check(sim_labels_counts, real_labels_counts, paste0(sim_name," counts"), paste0(real_name,"  counts"))
  if (!err_len$status) stop(err_len$msg)
  real_labels_counts<- real_labels_counts[sort(names(real_labels_counts))]
  sim_labels_counts<- sim_labels_counts[sort(names(sim_labels_counts))]
  OK <- ( sum(real_labels_counts != sim_labels_counts) ==0)
  msg <- NULL
  if (!OK){
    if (is.null(real_name)) msg<- paste0("Error: left and right have different histograms")
    else msg<-paste0("Error: ",sim_name,"  and ", real_name,"  have different histograms")
  }
  list(status=OK, msg=msg )
}

pick_matching_idx<- function(lab, real_labels){
  which(real_labels==lab)
}

reorder_left_by_right<- function(sim_labels, real_labels){
  "returns permutation of left argument that would reflect ordering of right argument,
  where it is assumed that left and right have same length, same members but differ only
  in order. "
  err_check <- len_check(sim_labels, real_labels, "simulated data labels","real data labels")
  if (!err_check$status) stop(err_check$msg)
  same_mem_check <- same_members_check(
    unique(sim_labels), unique(real_labels), "simulated data labels", "real data labels")
  if(!same_mem_check$status) stop(same_mem_check$msg)
  mem_c_check <-  members_count_check(sim_labels, real_labels,"simulated data labels", "real data labels")
  if(!mem_c_check$status) stop(mem_c_check$msg)
  distinct_labels <- unique(real_labels)
  label_idx_map <- lapply(distinct_labels, pick_matching_idx, real_labels)
  names(label_idx_map)<-distinct_labels
  result_perm <- c(1: length(sim_labels))
  for (distinct_lab in names(label_idx_map)){
   result_perm[label_idx_map[[distinct_lab]]] =  which(sim_labels==distinct_lab)
  }
  return(result_perm)
}




ReorderSimByReal<- function(sim_result, real_clustering_result ){
  #change order of elements in  simulated data
  #to reflect the order of elements in real data
  #this will work only if simulated data modules share names and gene counts with
  #modules of real data
  #reordering simulated data simplifies comparsions of statistics of simulated and real data,
  #as well as comparsions of clustering of the simulated data vs clustering of real data
  #sim_result is output from simulateDatExpr
  #real_clustering_result is in the form of output of ClusteringResults functions 
  #
  #function returns adjusted sim_result, that reflects ordering created in clustering real data
  real_labels<-real_clustering_result$color_labels
  sim_labels<-labels2colors(sim_result$allLabels)
  new_sim_order <- reorder_left_by_right(sim_labels, real_labels)
  adjusted_sdE             <- sim_result$datExpr[,new_sim_order]
  adjusted_setLabels       <- sim_result$setLabels[new_sim_order]
  adjusted_allLabels       <- sim_result$allLabels[new_sim_order]
  adjusted_trueKME         <- sim_result$trueKME[new_sim_order]
  adjusted_trueKME.whichMod<- sim_result$trueKME.whichMod[new_sim_order]
  list(datExpr=adjusted_sdE,
       setLabels=adjusted_setLabels,
       allLabels=adjusted_allLabels,
       labelOrder=sim_result$labelOrder,
       trueKME=adjusted_trueKME,
       trueKME.whichMod=adjusted_trueKME.whichMod)
}

BaseClustering2SavedSim<- function(new_sim_name,
                           datasets_path, base_dataset_name, base_network_name, base_clustering_name, 
                            ReorderByReal=TRUE,
                            base_expr_data_path=NULL, base_expr_data=NULL, #one of those must not be null
                           base_has_decision=FALSE, base_expr_RData=TRUE,
                           method_cor='spearman', #for adjaceny matrix
                           verbose=0, other_named_args=NULL,#simulateDatExpr_fromInput
                           min_from_perm=FALSE, abs_cor=TRUE, true_grey_frac=0.5, #clust_res2simDEargs
 ,...#clust_res2simDEargs
    ){
    clres<- DoClusterFromFilenameArgs(datasets_path, base_dataset_name, base_network_name, base_clustering_name, method_cor, base_expr_data, base_expr_data_path, base_has_decision, base_expr_RData, calculateMEs=TRUE)
if (is.null(base_expr_data)){ if (is.null(base_expr_data_path)) stop('expr_data or expr_data_path must not be null') 
                                                      else if (base_expr_RData) { load(base_expr_data_path); base_expr_data<-as.matrix(data.train) } 
                                                                            else  base_expr_data<- readRDS(base_expr_data_path)
                        }
if (base_has_decision) base_expr_data<- base_expr_data[,-1]
C<- cor(base_expr_data, method=method_cor)
  input_args_list= clust_res2simDEargs(base_expr_data, clres$color_labels, C, clres$ME_list$MEs, min_from_perm, abs_cor, true_grey_frac, save_steps=FALSE, save_final=FALSE, ...) 
 simul<-simulateDatExpr_fromInput(input_args_list, verbose, other_named_args)
 if (ReorderByReal)
    s_result<-ReorderSimByReal(simul$sim_result, clres) 
 else
    s_result<- simul$sim_result
    specs<- list(name=new_sim_name,
                 base_dataset_name=base_dataset_name,
                 base_network_name=base_network_name,
                 base_clustering_name=base_clustering_name,
                 reordered=ReorderByReal,
                 min_from_perm=min_from_perm,
                 abs_cor=abs_cor,
                 true_grey_frac=true_grey_frac)
    if (!is.null(other_named_args)) specs$other_named_args=other_named_args 
    prefix_path=paste0(datasets_path,'/',base_dataset_name,'/')
    dir.create(paste0(prefix_path,'simulations/base/',new_sim_name), showWarnings=FALSE, recursive=TRUE)
    saveRDS(specs, paste0(prefix_path,'simulations/base/',new_sim_name,'/Specs.rds'))
    return(s_result)
}

Specs2Sim<- function(specs,datasets_path,
                           base_expr_data_path=NULL, base_expr_data=NULL, #one of those must not be null
                           base_has_decision=FALSE, base_expr_RData=TRUE,
                           method_cor='spearman',verbose=0,...)
{
    clres<- DoClusterFromFilenameArgs(datasets_path, specs$base_dataset_name, specs$base_network_name, specs$base_clustering_name, method_cor, base_expr_data, base_expr_data_path, base_has_decision, base_expr_RData, calculateMEs=TRUE)
if (is.null(base_expr_data)){ if (is.null(base_expr_data_path)) stop('expr_data or expr_data_path must not be null') 
                            else if (base_expr_RData) { load(base_expr_data_path); base_expr_data<-as.matrix(data.train) }                              else  base_expr_data<- readRDS(base_expr_data_path)
                            }
if (base_has_decision) base_expr_data<- base_expr_data[,-1]
C<- cor(base_expr_data, method=method_cor)
  input_args_list= clust_res2simDEargs(base_expr_data, clres$color_labels, C, clres$ME_list$MEs, specs$min_from_perm, specs$abs_cor, specs$true_grey_frac, save_steps=FALSE, save_final=FALSE, ...) 
 if (is.null(specs$other_named_args))
 simul<-simulateDatExpr_fromInput(input_args_list, verbose)
 else
 simul<-simulateDatExpr_fromInput(input_args_list, verbose, other_named_args)
 if (specs$reordered)
    s_result<-ReorderSimByReal(simul$sim_result, clres) 
 else
    s_result<- simul$sim_result
 return(s_result)
}

ReadSimSpecsFile <- function(sim_name, datasets_path,dataset_name)
{
    prefix_path=paste0(datasets_path,'/',dataset_name,'/')
    MySpecs<-readRDS(paste0(prefix_path,'simulations/base/',sim_name,'/Specs.rds'))
    return(MySpecs)
} 

SpecsFile2Sim<- function(sim_name,datasets_path,dataset_name,
                           base_expr_data_path=NULL, base_expr_data=NULL, #one of those must not be null
                           base_has_decision=FALSE, base_expr_RData=TRUE,
                           method_cor='spearman',verbose=0,...)
{
    specs<- ReadSimSpecsFile(sim_name, datasets_path,dataset_name)
    return(Specs2Sim(specs,datasets_path,
                     base_expr_data_path, base_expr_data,
                           base_has_decision, base_expr_RData,
                           method_cor,verbose,...))
}
fromInput_simulate2ClusteringResults<- function(input_args_list,other_named_args, real_clustering_result,
                                                method_cor="spearman", adj_power=5, 
                                                save_steps=FALSE, save_final=TRUE,
                                                diss_as_TOM= TRUE,ts_other_args=NULL,
                                                hclust_method = "average", minModSize=30,
                                                dc_sensivity= 2, dc_pam2parent=FALSE, dc_doPAM = TRUE, 
                                                dc_other_args=NULL, ME_other_args= NULL,
                                                verbose=0) {
  #does simulateDatExpr_fromInput -> ReorderSimByReal -> cor -> adj -> ClusteringResults.fromAdjacency
  #objects from intermediate steps are not returned (excluding simulated expression data)
  #to spare memory, but are saved on hard drive if save_steps=TRUE 
  simul<-simulateDatExpr_fromInput(input_args_list, verbose=verbose, other_named_args)
  s_result<-ReorderSimByReal(simul$sim_result, real_clustering_result)
  rm(simul)
  s_cor<- cor(s_result$datExpr, method=method_cor)
  s_adj<- adjacency.fromSimilarity(s_cor, power=adj_power)
  rm(s_cor)
  s_clres<- ClusteringResults.fromAdjacency(s_adj, s_result$datExpr, save_steps, save_final,
              diss_as_TOM,ts_other_args,
              hclust_method, minModSize,
              dc_sensivity, dc_pam2parent, dc_doPAM, 
              verbose_lvl=verbose, dc_other_args, ME_other_args)
s_expr_data<-s_result$datExpr
rm(s_adj); rm(s_result)
list(datExpr=s_expr_data,
     clustering_result=s_clres)
}

