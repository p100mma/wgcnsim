FROM rocker/r-ver:4.0.3

COPY *.R root/

RUN install2.r --error --skipinstalled --ncpus -1 \
    matrixStats \
    Hmisc \
    splines \
	 foreach \
	 doParallel \
	 fastcluster \
 	 dynamicTreeCut \
 	 survival \
 	 BiocManager 


RUN Rscript -e "BiocManager::install('Biobase');"
RUN Rscript -e "BiocManager::install('GO.db');"
RUN Rscript -e "BiocManager::install('preprocessCore');"
RUN Rscript -e "BiocManager::install('impute');"
RUN Rscript -e "BiocManager::install('WGCNA');"
CMD ["R"]
