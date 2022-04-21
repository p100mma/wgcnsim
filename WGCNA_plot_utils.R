
standard_args_PlotDendroAndColors<- function(clust_result, side= "DTC",main="Gene dendogram and modules", text_colors=TRUE, ...){
  text_arg<-NULL
  if (text_colors) text_arg<- clust_result$color_labels
  plotDendroAndColors(clust_result$geneTree, clust_result$color_labels, side, dendroLabels=FALSE,
                      hang=0.03, rowText=text_arg, addGuide = TRUE, guideHang = 0.05, main=main, ... )
}