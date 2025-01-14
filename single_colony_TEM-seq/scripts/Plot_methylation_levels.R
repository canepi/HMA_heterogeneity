##########################################################################
# Purpose: R Script to plot 'Mean Methylation level' from metadata (.tsv)
# Output: Two png files displaying Mean methylation by batch and replicate.
#
# Date: 25.Jan.22
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
# To do:
##########################################################################

## Using metadata tsv which contains passing QC info.

plot_MeanMeth <- function(meta.data.tsv = "./2_results/colony/MOLM/cell_metadata.tsv",
                          output.dir ="figures/colony/MOLM",
                          Methylation.cut.off = 70,
                          Passing.QC.only = "TRUE" # If 'FASE' this will plot all samples.
                          ){
  
  
  
  #####################
  mkdirs <- function(fp) {
    if(!file.exists(fp)) {
      mkdirs(dirname(fp))
      dir.create(fp)
    }
  } 
  
  mkdirs(output.dir)
  ###################
  library(dplyr)
  library(readr)
  library(ggplot2)
  
  ###
  # Input variables
  fig.dir <- output.dir
  
  Theshold <- Methylation.cut.off
  ##################
  # Read in data
  colony_info <- read.csv(meta.data.tsv)
  
  if(Passing.QC.only == "TRUE"){
    colony_info <- colony_info %>% filter(Passed_seqmonk_QC %in% c("y","Y"))
  }
  
  colony_info %>% group_by(Batch, Treatment)
  
  # Plot 1
  meth_plot_by_batch <-   ggplot(colony_info, aes(x=Treatment, y = meanMeth, fill = Batch) ) +
    geom_violin(position=position_dodge(0.7), drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 100 , 10), limits = c(0,100)) + 
    geom_dotplot(binaxis='y', stackdir='center',
                 position=position_dodge(0.7), binwidth = 0.8, drop = FALSE) + 
    theme_minimal() +
    geom_hline(yintercept = Theshold,
               color = "darkred",
               linetype = "dashed", lwd = 0.5) +
    ggplot2::annotate("text", 
                      "unt", 
                      Theshold + 1.5, 
                      label = "Methylation threshold",
                      color = "darkred", size = 3)
  
  # Plot 2 (combined batches)
  
  meth_plot_combined_batch <-   ggplot(colony_info, aes(x=Treatment, y = meanMeth) ) +
    geom_violin(position=position_dodge(0.7)) +
    scale_y_continuous(breaks = seq(0, 100 , 10), limits = c(0,100)) + 
    geom_dotplot(binaxis='y', stackdir='center',
                 position=position_dodge(0.7), binwidth = 0.6, alpha = 0.4) + 
    theme_minimal() +
    geom_hline(yintercept = Theshold,
               color = "darkred",
               linetype = "dashed", lwd = 0.5) +
    ggplot2::annotate("text", 
                      "unt", 
                      Theshold + 1.5, 
                      label = "Methylation threshold",
                      color = "darkred", size = 3)
  
  # Save plots
  ggsave(meth_plot_by_batch, filename = paste0(fig.dir, "/Colony_methylation_levels_by_batches_and_Tx.png"))
  ggsave(meth_plot_combined_batch, filename = paste0(fig.dir, "/Colony_methylation_levels_combined_batches.png"))
  
}
