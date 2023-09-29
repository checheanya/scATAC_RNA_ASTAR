#install.packages('dplyr')
#install.packages('plyr')
#install.packages('ggplot2')
#install.packages('patchwork')

# run from terminal:  sudo apt-get install libcurl4-openssl-dev libxml2-dev

#install.packages("RCurl")
#install.packages("XML")

#install.packages("Seurat")
#install.packages("Signac")
#install.packages('cowplot')
#install.packages('hdf5r')

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.17")
#BiocManager::install("EnsDb.Hsapiens.v86")
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#BiocManager::install("biovizBase")
#BiocManager::install("glmGamPoi")
#BiocManager::install("GenomicRanges")

# something in downgrading Signac to it throws layer error...
remotes::install_github("stuart-lab/signac", "seurat5")

# do:
# sudo add-apt-repository ppa:ubuntugis/ppa
# sudo apt-get update
# sudo apt -y install libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev  
# sudo apt -y install libgsl-dev
# sudo apt-get install build-essential r-cran-raster libgdal-dev gdal-bin

# install.packages("devtools")

# if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
#library(devtools)
#install.packages('doSNOW')
#install.packages('plot3D')
#install_github("r3fang/SnapATAC")

#devtools::install_github("aertslab/RcisTarget")
#devtools::install_github("aertslab/AUCell")
#devtools::install_github("aertslab/cisTopic")


# STATISTICS
#install.packages('entropy')
#install.packages('lme4')

# ARCH R
#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
#library(ArchR)
#ArchR::installExtraPackages()




















