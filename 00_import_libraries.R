################################################################################
# This script is used to clean the envrionment and import all necessary packages
################################################################################
# Specify the list of packages that need to be imported ########################
list.of.packages <- c("Seurat",
                      "SingleCellExperiment",
                      "optparse", 
                      "comprehenr", 
                      "tidyverse", 
                      "ggplot2", 
                      "SoupX",
                      "comprehenr",
                      "DoubletFinder",
                      "vroom",
                      "hash",
                      "DT",
                      "janitor",
                      "knitr",
                      "circlize",
                      "formattable",
                      "htmlwidgets",
                      "plotly",
                      "stringr"
)

bioc.packages <- c("celda", 
                   "BiocSingular", 
                   "PCAtools", 
                   "SingleCellExperiment",
                   "scRepertoire", 
                   "sctransform", 
                   "org.Mm.eg.db")

# Check if packages are installed ##############################################

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.bioc.packages <- bioc.packages[!(bioc.packages %in% installed.packages()[,"Package"])]

# Install new packages #########################################################
if(length(new.packages)) install.packages(new.packages)
if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages, update = FALSE)

# Import all packages ##########################################################
package_loading_Status <- lapply(list.of.packages, 
                                 require, 
                                 character.only = TRUE)

package_loading_Status_bioc <- lapply(bioc.packages, 
                                      require, 
                                      character.only = TRUE)

##### SPECIAL UPDATE FOR THE PACKAGE CLUSTER PROFILER
if ("clusterProfiler" %in% installed.packages()[, "Package"] == FALSE) {
  devtools::install_github("YuLab-SMU/clusterProfiler", upgrade = "never")
} else {
  if (packageVersion("clusterProfiler") != "4.7.1.3" ){
    remove.packages("DOSE")
    remove.packages("clusterProfiler")
    devtools::install_github("YuLab-SMU/clusterProfiler", upgrade = "never")
  }
}
library("clusterProfiler")

# EOF ##########################################################################
