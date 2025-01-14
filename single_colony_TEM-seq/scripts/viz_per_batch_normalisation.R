##########################################################################
# Purpose: R Script to visualise 'per_batch_nornmalisation.R'
# Output: Produces PCA plots for...
#
# Date: 25.August.22
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: This was originally developed as part of a package build running on unix/HPC by Al, which would allow shell style submission.
# To do: Embed directly within per_batch_normalisation.R script
##########################################################################
## Load packages
source("scripts/converted_from_HPC/library-viz.R")
source("./scripts/converted_from_HPC/plot_reducedDim.R") # Custom plotting function loaded.

viz_per_batch_normalisation <- function(output.file = '2_results/colony/HL60/Batch_Correction/sce_qc2_PerBatchNormPCA.rds', # output RDS path to save the data
                                        output.fig = 'figures/colony/HL60/Batch_Correction', # output folder for figures
                                        output.fig.name = "PerBatch_QC2_PCA.png", # output name for figures
                                        shape = 1,
                                        size = 3,
                                        width = 9,
                                        height = 4,
                                        all_cols = c("unt" = "green", "aza" =  "blue", "dac"  = "red"), # Change as required for your batch c("rep_name1" = "colour1, "rep_name2 = "colour2")
                                        reducedDim = 'PCA.batch'){

  ## viz per-batch normalisation

  cat("\nloading the output file ...\n")
  message(output.file)
  sce.batch <- readRDS(output.file)
  
  plot_reducedDim(
    sce = sce.batch,
    reducedDim = reducedDim,
    comps = c(1, 2),
    colBy = 'Treatment', facet_by = 'Batch',
    cols = all_cols,
    axis.title = 'PC',
    size = size,
    shape = shape)

  filename <- here(paste0(output.fig, "/", output.fig.name))
  cat("\nwriting: ", filename, "\n")
  ggsave(filename =  filename,
         width = width,
         height = height)
}
