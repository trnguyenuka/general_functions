gc()
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
wd <- getwd()
source(file.path(wd, "00_import_libraries.R"))
source(file.path(wd, "00_helper_functions.R"))

#####----------------------------------------------------------------------#####
# INPUT ARGUMENTS
#####----------------------------------------------------------------------#####
path.to.input <- file.path(wd, "input")

path.to.output <- file.path(wd, "output")
dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)

path.to.html <- file.path(wd, "html")
dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)

path.to.template <- file.path(wd, "templates")

clusters.to.be.removed <- c(16)

name.of.s.obj <- "1st_dataset_removed_7_9.without_reInt.rds" 

with.re.integration <- FALSE

if (with.re.integration == TRUE){
  status <- "with_reInt"
} else if (with.re.integration == FALSE){
  status <- "without_reInt"
} else {
  stop("The parameter with.re.integration must be either TRUE or FALSE")
}

save.dataset.name <- "1st_dataset_removed_7_9"  # <<<<< CHANGE HERE

save.obj.name <- sprintf("%s_removed_%s.%s.rds", save.dataset.name, paste(clusters.to.be.removed, collapse = "_"), status)
save.html.name <- str_replace(save.obj.name, ".rds", ".html")

path.to.output <- file.path(path.to.output, str_replace(save.obj.name, ".rds", ""))
dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)

##### QUICK CHECK
fetch.available.objs <- basename(Sys.glob(file.path(path.to.input, "*")))

if (with.re.integration %in% c(TRUE, FALSE) == FALSE){
  stop("The parameter with.re.integration must be TRUE or FALSE.")
}

if (name.of.s.obj %in% fetch.available.objs == FALSE){
  stop("Specified dataset not available in the input folder. Please choose another dataset")
}

if (is.null(clusters.to.be.removed) == TRUE){
  stop("Input empty! Please insert a list of clusters to be removed")
}

# Fixed parameters
chosen.seed <- 42
num.dim.integration <- 25 
num.PCA <- 25
num.dim.cluster <- 25
cluster.resolution <- 1
num.PC.used.in.Clustering <- 25

#####----------------------------------------------------------------------#####
# MAIN SCRIPTS
#####----------------------------------------------------------------------#####

##### 1. Processing the Seurat object
s.obj <- readRDS(file.path(path.to.input, name.of.s.obj))

s.obj.removed <- subset(s.obj, seurat_clusters %in% clusters.to.be.removed == FALSE)
DefaultAssay(s.obj.removed) <- "RNA"

if (with.re.integration == TRUE){
  
  data.list <- SplitObject(s.obj.removed, split.by = "name")
  data.list <- lapply(X = data.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})
  
  k.filter <- 200
  
  anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:num.dim.integration, scale=F,
                                    k.filter = k.filter) ## THIS IS CCA DIMENSIONS
  
  s.obj_inte <- IntegrateData(anchorset = anchors, dims = 1:num.dim.integration, k.weight = k.filter) ## THIS IS PCA DIMENSION
  
  s.obj_inte <- s.obj_inte[, colnames(s.obj.removed)]
  
  s.obj.removed[['integrated']] <- s.obj_inte[['integrated']]
  
  s.obj.removed@commands <- c(s.obj.removed@commands, s.obj_inte@commands)
  
  s.obj.removed@tools <- c(s.obj.removed@tools, s.obj_inte@tools)
  
  DefaultAssay(s.obj.removed) <- "integrated"
  
  remove.YFP.VariableFeatures <- to_vec(for (item in VariableFeatures(s.obj.removed)) if (item != "YFP") item)
  
  s.obj.removed <- ScaleData(s.obj.removed, verbose = FALSE, features = remove.YFP.VariableFeatures)
  
  s.obj.removed <- RunPCA(s.obj.removed, npcs = num.PCA, verbose = FALSE, reduction.name="INTE_PCA")
  
  s.obj.removed <- RunUMAP(s.obj.removed, reduction = "INTE_PCA", dims = 1:num.PCA, reduction.name="INTE_UMAP")
  
  s.obj.removed <- FindNeighbors(s.obj.removed, reduction = "INTE_PCA", dims = 1:num.dim.cluster)
  
  s.obj.removed <- FindClusters(s.obj.removed, resolution = cluster.resolution)
  
} else {
  status <- "without_reInt"
  
  s.obj.removed <- NormalizeData(s.obj.removed) # ---> use Log Normalized
  s.obj.removed <- FindVariableFeatures(s.obj.removed, selection.method = "vst")
  
  remove.YFP.VariableFeatures <- to_vec(for (item in VariableFeatures(s.obj.removed)) if (item != "YFP") item)
  
  s.obj.removed <- ScaleData(s.obj.removed, features = to_vec(for (item in rownames(s.obj.removed)) if(item != "YFP") item))
  
  s.obj.removed <- RunPCA(s.obj.removed, npcs = num.PCA, verbose = FALSE, reduction.name="INTE_PCA", features = remove.YFP.VariableFeatures)
  
  s.obj.removed <- RunUMAP(s.obj.removed, reduction = "INTE_PCA", dims = 1:num.PCA, reduction.name="INTE_UMAP", seed.use = chosen.seed)
  # s.obj.removed <- RunTSNE(s.obj.removed, reduction = "INTE_PCA", dims = 1:num.PCA, reduction.name="INTE_TSNE", seed.use = chosen.seed)
  
  s.obj.removed <- FindNeighbors(s.obj.removed, reduction = "INTE_PCA", dims = 1:num.PC.used.in.Clustering)
  s.obj.removed <- FindClusters(s.obj.removed, resolution = cluster.resolution, random.seed = chosen.seed)
}

saveRDS(s.obj.removed, file.path(path.to.output, save.obj.name))

##### 2. Generate report from defined templates
rmarkdown::render(input = file.path(path.to.template, "01_generate_cluster_DE_genes.Rmd"), 
                  params = list(
                    path.to.input.sobj = file.path(path.to.output, save.obj.name),
                    path.to.wd = wd,
                    dataset_name = save.obj.name,
                    path.to.output = path.to.output
                  ),
                  output_file = save.html.name,
                  output_dir = path.to.html)  



