##########################################################################
# Purpose: R Script to visualise variance explanation from sce.RDS files.
# Output: PCA plot of MNN corrected data
#
# Date: 25.August.22
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: This was originally developed as part of a package build running on unix/HPC by Al, which would allow shell style submission.
# To do: Update plot size depending on number of input files, or enable as an input variable.
##########################################################################
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(ggplot2, knitr, meta, scuttle, cowplot, dplyr, scater, here)



#
viz_variance_explained <- function(sce.rds = c("./2_results/colony/HL60/Batch_Correction/sce_qc1_PerBatchNormPCA.rds", 
                                               "./2_results/colony/HL60/Batch_Correction/sce_qc2_PerBatchNormPCA.rds"),
                                  exprs.values.name = c("counts","counts"), # Needs to be the same length as the number of sce.rds files.
                                  plot.names = c("Raw data (QC1)", "QC filtered (QC2)"),
                                  output.file.dir = "figures/colony/HL60/Batch_Correction",
                                  output.file.name = "colony-explanatoryVariables-mnncorrected.png",
                                  variables.to.measure = c("total", "Batch", "Treatment", "subsets_Gene_detected", "subsets_TE_detected", "subsets_Mito_detected")){

  cat("Reading in",length(sce.rds), "rds files: \n", 
      sce.rds)
  cat("\nAnd extracting/plotting variance explained by: ", variables.to.measure)
  
  i <- 0
  for(rds in sce.rds){
    i <- i+1
    
    rds.data <- readRDS(rds)
    
    varianace_explained <- getVarianceExplained(rds.data, exprs_values = exprs.values.name[i], 
                                                variables = variables.to.measure)
    
    assign(paste0("plot_", i),
           plotExplanatoryVariables(varianace_explained) +
             ggtitle(plot.names[i]) +
             theme(plot.title = element_text(hjust = 0.5)))
    
  }
  
  
  plot_to_save <- plot_1 # This will be overwritten if more than 1 file has been provided (see 'if' function below).
  
  if(length(sce.rds) == 2){
    plot_to_save <- 
      plot_grid(nrow = 2,
                plot_1,
                plot_2)
  } else if(length(sce.rds) == 3) {
    plot_to_save <- 
      plot_grid(nrow = 3,
                plot_1,
                plot_2,
                plot_3)
  } else if(length(sce.rds) == 4) {
    plot_to_save <- 
      plot_grid(nrow = 4,
                plot_1,
                plot_2,
                plot_3,
                plot_4)
  } else if(length(sce.rds) == 5) {
    plot_to_save <- 
      plot_grid(nrow = 5,
                plot_1,
                plot_2,
                plot_3,
                plot_4,
                plot_5)
  } else if(length(sce.rds) == 6) {
    plot_to_save <- 
      plot_grid(nrow = 6,
                plot_1,
                plot_2,
                plot_3,
                plot_4,
                plot_5,
                plot_6)
  }
  
  cat("\n",length(sce.rds), "plots being saved to:\n")
  message(paste0(output.file.dir, "/",output.file.name))
  
  ggsave(plot_to_save, 
         filename = paste0(output.file.dir, "/",output.file.name), width = 7, height = 12, bg = "white")
  
}
################################################################