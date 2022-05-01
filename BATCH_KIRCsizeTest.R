source('initial.R')
tm<-readRDS("KIRCTOM.rds")
dg<-readRDS("KIRCdg.rds")

sizes<- floor( rep( 20221, 9) / 2^9)
labelings<- vector(mode="list", length=9)
for (i in seq_along(  sizes )){
    print(i)
   labelings[[i]]<-labels2colors(dendro2clusters(dg, tm, doPAM=FALSE, minModuleSize=sizes[[i]]) i)
}

saveRDS(labelings, 'KIRChierMAJOR.rds')

rm(tm)
for (i in seq_along(  sizes )){
    print(i)
   jpeg(paste0("KIRCsize",size,".jpg"),width=1000, height=600)
   standard_args_PlotDendroAndColors(list(geneTree=dg,
                            color_labels=labelings[[i]]), main=paste0("minSize",size), side=size)
    dev.off()
}


   jpeg(paste0("KIRCALL",".jpg"),width=1000, height=900)
   plotDendroAndColors(dg, as.data.frame(labelings), 
   sizes)                    
    dev.off()

