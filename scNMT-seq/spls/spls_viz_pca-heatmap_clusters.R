##########################################################################
# Purpose: R Script to visualise PCAs from MixOmics (s)PLS functions.
# Output: PCA plots (png)
#
# Date: 15.08.2023
# Version: v.0.1.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) 
# To do: 
## 1) 
# Websites:
## 
##########################################################################
# Load packages
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(ggplot2, dplyr, knitr, SingleCellExperiment, rlang, Seurat, gsubfn, data.table, cowplot, gridExtra, grid, RColorBrewer, circlize)
source("./scripts/utils_custom.R")


#########################################################################################################################################################
# Variables
# Input files
input.spls.file = "./2_results/scNMT/NMT/Autosomal/pct10_Y100-100_Tx_only/mint.block.spls_NA_rna_rna_pctCells_10_nhvrs_0_ncomps_3_mode_canonical.rds"
#cell.metadata.file = "./2_results/scNMT/NMT/Autosomal/cell_metadata_updated.tsv"
clusters.file = "./figures/scNMT/NMT/Autosomal/pct10_Y100-100_Tx_only/heatmap/Heatmap_NMT_C1-2_RNA_rowKmeans3_colKmeans4-ClusterIDs_SAMPLES.csv"


# Output
output.fig.dir = "./figures/scNMT/NMT/Autosomal/pct10_Y100-100_Tx_only/heatmap/spls_comps/"

# Additional variables
#treatment_colours = c(DAC = "#00FFFF", 
#                      #Unt = "orange", 
#                      AZA = "orchid")
plot.width = 10
plot.height = 10
block.to.plot =  "rna"
#########################################################################################################################################################

# Create dir
dir.create(output.fig.dir, showWarnings = FALSE, recursive = TRUE)

# Import data
## Read in input files
cat("\nReading in spls res file: ", input.spls.file)
res <- readRDS(input.spls.file)
#cat("\nReading in cell metadata file: ", cell.metadata.file)
#cell.metadata <- fread(cell.metadata.file) %>% as.data.frame()
cat("\nReading in clusters file: \n", clusters.file)
clusters <- fread(clusters.file, colClasses = c("character","factor","factor","factor")) 
# Ensure factors are in correct order
clusters$Batch <- factor(clusters$Batch, levels = c("b320", "b620", "b322"))
clusters$cluster <- factor(clusters$cluster, levels = c(1,2,3,4))



cat("\nextracting variates ...\n")
variates <- lapply(res$variates, function(x){
  x <- as.data.table(x, keep.rownames = 'sample')
})

if (!('rna' %in% names(variates)))
{
  ## in spls analysis name Y as 'rna'
  names(variates)[which(names(variates) == 'Y')] <- 'rna'
}

variates <- rbind_df_list(variates, new_col = 'assay')
variates <- as.data.table(variates, keep.rownames = 'sample')

plots_to_save <- NULL
all_blocks <- c( "met_CGI","met_H3K27ac", "met_H3K4me3", "met_Promoter", "met_Window3k1k", 
                 "acc_CGI", "acc_H3K27ac", "acc_H3K4me3", "acc_Promoter", "acc_Window3k1k",
                 "rna")
for(block in all_blocks){
  
  block.to.plot <- block
  # Filter for block of interest to plot
  cat("\nFiltering variates and plotting only the context: ", block.to.plot)
  variates_sub <- variates %>% 
    filter(assay == block.to.plot)
  
  # Combine with cluster table
  variates_sub <- variates_sub %>%
    left_join(clusters, by = c("sample" = "sample_ID")) 
  
  
  plots_to_save[[block]] <- 
    ggplot(variates_sub, aes(x = comp1, y = comp2, col = cluster)) +
    geom_point(aes(shape = Treatment), alpha = 0.7) + 
    #scale_colour_gradientn(colors = my_colors) +
    theme_classic() +  # + stat_ellipse(type = "norm", linetype = 2,level = 0.8) +
    ggtitle(block) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
}

# Save plots 
## Save rna as svg
block <- "rna"
output.fig.name <- paste0(output.fig.dir,block,"-heatmap_Kmeans4.svg")
cat("\n Saving: ", output.fig.name)
ggsave(plots_to_save[[block]], filename = output.fig.name, height = 5, width =6)

## Individual
for(block in all_blocks) {
  output.fig.name <- paste0(output.fig.dir,block,"-heatmap_Kmeans4.png")
  cat("\n Saving: ", output.fig.name)
  ggsave(plots_to_save[[block]], filename = output.fig.name, height = 5, width =6)
}

## combined (cow)plots
### Met
output.fig.name <- paste0(output.fig.dir,"Met_combined_plots-heatmap_Kmeans4.png")
cat("\n Saving combined plot: ", output.fig.name)
plot_grid(plots_to_save[[1]],plots_to_save[[2]], plots_to_save[[3]],
          plots_to_save[[4]],plots_to_save[[5]], plots_to_save[[11]],
          ncol = 3, nrow =2)
ggsave(filename = output.fig.name, height = 10, width = 18)
### Acc
output.fig.name <- paste0(output.fig.dir,"Acc_combined_plots-heatmap_Kmeans4.png")
cat("\n Saving combined plot: ", output.fig.name)
plot_grid(plots_to_save[[6]],plots_to_save[[7]], plots_to_save[[8]],
          plots_to_save[[9]],plots_to_save[[10]], plots_to_save[[11]],
          ncol = 3, nrow =2)
ggsave(filename = output.fig.name, height = 10, width = 18)

cat("Finished! For now....")

