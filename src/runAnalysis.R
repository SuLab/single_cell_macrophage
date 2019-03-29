#### Run Analysis ####

### paths must be set first
## repository directory (e.g. ~/something/something/single_cell_macrophage)
repo_dir <- "/single_cell_macrophage"

#### render to html ####
## set controller script and output paths
# ctrl_path <- paste(repo_dir, "/src/controller/macrophage_all.Rmd", sep = "")
# 
# date <- gsub("-", "_", Sys.Date())
# output_path <- paste(repo_dir, "/results/macrophage_analysis_", date, ".html", sep = "")
# 
# ## render
# rmarkdown::render(
#   input = ctrl_path, 
#   output_format = "html_document", 
#   output_file = output_path)


#### render to pdf ####
## set controller script and output paths
ctrl_path <- paste(repo_dir, "/src/controller/macrophage_all.Rmd", sep = "")

date <- gsub("-", "_", Sys.Date())
output_path <- paste(repo_dir, "/results/macrophage_analysis_", date, ".pdf", sep = "")

## render
rmarkdown::render(
  input = ctrl_path,
  output_format = "pdf_document",
  output_file = output_path)