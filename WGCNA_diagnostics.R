

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


### Simulated vs Real comparsion
### descriptions of simulated data are assumed to be obtained after...
### ...clustering of simulated data in the same way that real data was clustered, unless stated otherwise.
### they are not necessairly equal to the assumed descriptions in the data generating process



RealSimContingency <- function(true_labels, sim_labels) {
  #true_lables <- colors of genes in real data
  #sim_labels <- colors of genes in simulated data clustered
  #in principle can be used with any equal length string vectors.
  table(true_labels, sim_labels)
}

RealSimRandIndex <- function(true_labels, sim_labels, adjust=FALSE){
  #compare similarity of clustering of real data and simulated data
  #by RandIndex measure. in principle can be used with any equal length string vectors.
  randIndex(RealSimContingency(true_labels,sim_labels), adjust=adjust)
}

RealSimEigenCor<- function(real_MEs, sim_MEs, abs_cor=FALSE, signif_digits=NULL, method="pearson" ){
  #output: rows - real data, columns - simulated
  res<-cor(real_MEs, sim_MEs, method=method)
  rownames(res)<-paste0("REAL_",rownames(res))
  colnames(res)<-paste0("SIM_",colnames(res))
  if (abs_cor) res<- abs(res)
  if (!is.null(signif_digits)) res <- signif(res, signif_digits)
  return(res)
}