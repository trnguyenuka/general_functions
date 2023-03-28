gc()
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
wd <- getwd()
source("00_import_libraries.R")
source("00_helper_functions.R")

#####----------------------------------------------------------------------#####
# INPUT ARGUMENTS
#####----------------------------------------------------------------------#####
path.to.input <- "./input"

path.to.output <- "./output"
dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)
path.to.html <- "./html"
dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)

path.to.template <- "./templates"

name.of.s.obj <- "merged_all_first_dataset_BCR.integrated.rds" 

cluster1 <- 4 # <<<<< CHANGE HERE
cluster2 <- 10 # <<<<< CHANGE HERE

# save object name
save.dataset.name <- "test"  # <<<<< CHANGE HERE

save.html.name <- sprintf("%s_%s.DEA.html", save.dataset.name, paste(c(cluster1, cluster2), collapse = "_"))

##### QUICK CHECK
fetch.available.objs <- basename(Sys.glob(file.path(path.to.input, "*")))

if (name.of.s.obj %in% fetch.available.objs == FALSE){
  stop("Specified dataset not available in the input folder. Please choose another dataset")
}

if (is.null(cluster1) == TRUE){
  stop("Input empty! Please insert the first cluster")
}

if (is.null(cluster2) == TRUE){
  stop("Input empty! Please insert the first cluster")
}


s.obj <- readRDS(file.path(path.to.input, name.of.s.obj))

if (cluster1 %in% unique(s.obj$seurat_clusters) == FALSE){
  stop(sprintf("The selected cluster 1 = %s does not exists in the data object", cluster1))
}
 
if (cluster2 %in% unique(s.obj$seurat_clusters) == FALSE){
  stop(sprintf("The selected cluster 2 = %s does not exists in the data object", cluster2))
}

s.obj <- subset(s.obj, seurat_clusters %in% c(cluster1, cluster2))

rmarkdown::render(input = file.path(path.to.template, "02_compare_any_2_clusters.Rmd"), 
                  params = list(
                    s.obj = s.obj, 
                    cluster1 = cluster1,
                    cluster2 = cluster2,
                    name_2_clusters = sprintf("Cluster %s", paste(c(cluster1, cluster2), collapse = ", ")),
                    path.to.wd = wd
                  ),
                  output_file = save.html.name,
                  output_dir = path.to.html)  


