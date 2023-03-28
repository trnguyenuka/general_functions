gc()
rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
wd <- getwd()
source("00_import_libraries.R")
source("00_helper_functions.R")
library("clusterProfiler")
#####----------------------------------------------------------------------#####
# INPUT ARGUMENTS
#####----------------------------------------------------------------------#####
path.to.input <- file.path(wd, "input")

path.to.output <- file.path(wd, "output")
dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)

path.to.html <- file.path(wd, "html")
dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
# input.gene.list.file <- "DE analysis between 2 clusters, Cluster 5, 10.csv" # <<<< CHANGE HERE
# input.gene.list.file <- "DE analysis between 2 clusters, Cluster 1, 2.xlsx"
# input.gene.list.file <- "test_cluster_4_vs_cluster_9.raw.xlsx"
input.gene.list.file <- "test_cluster_4_vs_cluster_10.raw.xlsx"
#####----------------------------------------------------------------------#####

if (grepl(".csv", input.gene.list.file) == TRUE){
  file.ext <- ".csv"
  inputdf <- read.csv(file.path(path.to.input, input.gene.list.file)) %>%
    subset(select = -c(X)) %>% arrange(desc(avg_log2FC))  
} else if (grepl(".xlsx", input.gene.list.file) == TRUE){
  file.ext <- ".xlsx"
  inputdf <- readxl::read_excel(file.path(path.to.input, input.gene.list.file))
  if ("DE analysis between" %in% colnames(inputdf)[[1]]) {
    inputdf <- readxl::read_excel(file.path(path.to.input, input.gene.list.file), skip = 1) %>%
      subset(select = -c(`...1`))
  }
}

path.to.template <- "./templates"

##### CHECK
if ("Gene" %in% colnames(inputdf) == FALSE | "avg_log2FC" %in% colnames(inputdf) == FALSE){
  stop("Gene or avg_log2FC column does not exists in the data input.")
}

inputdf <- inputdf %>% arrange(desc(avg_log2FC))

save.html.name <- sprintf("Pathway analysis %s", str_replace(input.gene.list.file, file.ext, ".html"))

rmarkdown::render(input = file.path(path.to.template, "03_pathway_analysis.Rmd"), 
                  params = list(
                    inputdf = inputdf
                  ),
                  output_file = save.html.name,
                  output_dir = path.to.html)  



