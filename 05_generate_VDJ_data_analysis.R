##### THIS SCRIPT APPLYIES TO THE [B.SIMONS] 1ST DATASET ONLY!
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

name.of.s.obj <- "1st_dataset_removed_7_9_and_16_and_9_removed_9.without_reInt.rds" 

# save object name
save.dataset.name <- "VDJ_removed_7_9_and_16_and_9"  # <<<<< CHANGE HERE

save.html.name <- sprintf("%s.VDJ.html", save.dataset.name)

##### QUICK CHECK
fetch.available.objs <- basename(Sys.glob(file.path(path.to.input, "*")))

if (name.of.s.obj %in% fetch.available.objs == FALSE){
  stop("Specified dataset not available in the input folder. Please choose another dataset")
}

s.obj <- readRDS(file.path(path.to.input, name.of.s.obj))

rmarkdown::render(input = file.path(path.to.template, "05_VDJ_data_1st_dataset.Rmd"), 
                  params = list(
                    s.obj = s.obj, 
                    path.to.wd = wd
                  ),
                  output_file = save.html.name,
                  output_dir = path.to.html)  


