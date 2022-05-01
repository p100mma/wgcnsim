source('initial.R')
tm<-readRDS("KIRCTOM.rds")
dg<-readRDS("KIRCdg.rds")

sizes<-seq(1200,3,-100) 
labelings<- vector(mode="list", length=length ( sizes ))
for (i in seq_along(  sizes )){
    print(i)
   labelings[[i]]<-labels2colors(dendro2clusters(dg, tm, doPAM=FALSE, minModuleSize=sizes[[i]]) )
}

saveRDS(labelings, 'KIRChierMID.rds')

rm(tm)
for (i in seq_along(  sizes )){
    print(i)
   jpeg(paste0("KIRCsize",sizes[[i]],".jpg"),width=1000, height=600)
   standard_args_PlotDendroAndColors(list(geneTree=dg,
                            color_labels=labelings[[i]]), main=paste0("minSize",sizes[[i]]), side=sizes[[i]])
    dev.off()
}


   jpeg(paste0("KIRCMID",".jpg"),width=1000, height=900)
   plotDendroAndColors(dg, as.data.frame(labelings), 
   sizes, dendroLabels=FALSE)                    
    dev.off()


tm<-readRDS("KIRCTOM.rds")
sizes<-seq(50,3,-7)
labelings<- vector(mode="list", length=length ( sizes ))
for (i in seq_along(  sizes )){
    print(i)
   labelings[[i]]<-labels2colors(dendro2clusters(dg, tm, doPAM=FALSE, minModuleSize=sizes[[i]]) )
}

saveRDS(labelings, 'KIRChierMINOR.rds')

rm(tm)
for (i in seq_along(  sizes )){
    print(i)
   jpeg(paste0("KIRCsize",sizes[[i]],".jpg"),width=1000, height=600)
   standard_args_PlotDendroAndColors(list(geneTree=dg,
                            color_labels=labelings[[i]]), main=paste0("minSize",sizes[[i]]), side=sizes[[i]])
    dev.off()
}


   jpeg(paste0("KIRCMINOR",".jpg"),width=1000, height=900)
   plotDendroAndColors(dg, as.data.frame(labelings), 
   sizes, dendroLabels=FALSE)                    
    dev.off()

