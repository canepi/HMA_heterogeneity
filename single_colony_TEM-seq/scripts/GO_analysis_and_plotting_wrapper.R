##########################################################################
# Purpose: R script wrapper go 'GO_analysis_and_plotting.R'
# Output: See 'GO_analysis_and_plotting.R'
#
# Date: 09.Feb.23
# Version: v.0.2.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## Added plot.key.colour which can be p.adjust, pvalue or q.value
## Added qvalue threshold options.
## Modified to handle inclusion of more than one cluster.number as input. 
## Added option for selecting multiple figure extension types
# To do: 
##########################################################################
source("./scripts/Go_analysis_and_plotting.R")


Perform_GO_and_plot_from_cluster_results <- function(Cluster.run = "./2_results/colony/HL60/DAC_high_vs_low/Methylation_Exp_mnnCorrected_correlation_and_heatmap_cluster_Run3.csv",
                                                     cluster.number = 1,
                                                     sce_used ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds", # For background
                                                     batch_correction_method = "mnnCorrected",
                                                     p.adj.method = "fdr",
                                                     p.adj.threshold = 0.05,
                                                     q.val.threshold = 0.4,
                                                     plot.key.colour = "p.adjust",
                                                     output.fig.dir = "./figures/colony/HL60/DAC_high_vs_low/",
                                                     output.file.dir = "./2_results/colony/HL60/DAC_high_vs_low/",
                                                     output.extension.name = "",
                                                     figure.extension = c("png", "svg")){ # Can be a single or multiple types. But beware. The figure sizes haven't been adjusted for the file type...


  # Need to manually extract ensembl hit list and background list as vectors. The function will then handle the rest!
  # Obtain hit list
  cat("\nFilter for genes in cluster ", cluster.number, " from '", Cluster.run, "'")
  cat("\nAnd filtering for genes")

  Cluster.run <-  read.csv(Cluster.run)
  ensembl_hit_list <- Cluster.run %>% 
    filter(cluster_number %in% cluster.number) %>% 
    dplyr::select(matches("genes|ensembl.TE")) %>% # This now flexibly allows selection of either column name. Could have just taken first column, but prefered to make sure the right column was selected (in case in a different order).
    .[[1]]
    # .[["ensembl.TE"]] 

  # Obtain background list
  sce_used <- readRDS(sce_used)
  sce_used_batch_correction <- assay(sce_used, batch_correction_method)
  background_list <- rownames(sce_used_batch_correction)



  Perform_Go_and_plot(ensembl.hit.list = ensembl_hit_list, # Needs to be a vector. 
                      background.list = background_list, # Needs to be a vector.
                      p.adj.method = p.adj.method, # Can be "NULL", "fdr","bh" etc
                      p.adj.threshold = p.adj.threshold, # Default 'always' being 0.05
                      q.val.threshold = q.val.threshold,
                      plot.key.colour = plot.key.colour, 
                      output.fig.dir = output.fig.dir,
                      output.file.dir = output.file.dir,
                      output.extension.name = output.extension.name, # This can be left blank with "", and all plots will remain with their basic names.
                      figure.extension = figure.extension)

}
