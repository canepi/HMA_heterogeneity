##########################################################################
# Purpose: R Script to plot comparison of RNAseq counts pre and post-correction.
# Output: ggplot of RNAseq counts averaged across features, split by treatment, batch and correction method (raw, normalised and corrected)
#
# Date: 14 Mar 2023
# Version: v.0.0.2
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 
# To do: 
##
# Weblinks:
##
#########################################################################
# Load packages
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(ggplot2, 
               dplyr,
               knitr, 
               stringr,
               tidyr,
               readr,
               data.table,
               SingleCellExperiment,
               rlang,
               scater,
               Seurat,
               gsubfn,
               cowplot)

plot_corrected_counts_comparison <- function(input.seurat.rds = "2_results/scNMT/RNAseq/Batch_Correction/Seurat/Seurat4.rds",
                                             output.fig = "figures/scNMT/RNAseq/Batch_Correction/Seurat/Seurat_correction_counts_comparison.png",
                                             batch.order = c("b320","b620","b322"),
                                             batch.names = c("060320" =  "b320", "250620" = "b620", "220331" = "b322"),
                                             raw.assay = "RNA",
                                             norm.assay = "SCT",
                                             corrected.assay = "integrated"){



  # Read in seurat object
  cat("Reading in: ", input.seurat.rds)
  seurat.corrected <- readRDS(input.seurat.rds)




  # RAW 
  cat("Reading, summarising and plotting feature counts for:", raw.assay, "assay")
  ## Convert counts to log2(x+1)
  Raw.df <-  GetAssayData(seurat.corrected[[raw.assay]]) %>%
  
    as.matrix(.) %>% + 1 %>% log2() %>% ## Convert counts to log2(x+1)
  
    reshape2::melt(., id = rownames(.)) %>% # Convert to long table for easier plotting
    dplyr::rename("id" = 1, "sample" = 2, "value" = 3) %>% 
    tidyr::separate(col ="sample", into = c('Well', 'Treatment', 'Batch'), remove = FALSE) %>% # This take a while due to so many columns created when melting each row to show genes by cell. If only after the average expression. Could probably average before...
    group_by(Batch, Treatment, id) %>% summarise(mean = mean(value))
  ## Update batch names and specify order for plotting
  Raw.df$Batch <- batch.names[Raw.df$Batch]
  Raw.df$Batch <- factor(Raw.df$Batch, level = batch.order)
  
  # plot

  Raw.df_plot <- 
    ggplot(data=Raw.df, aes(x=Batch, y=mean, fill=Batch)) + 
      geom_boxplot()  + 
      theme_bw() +
      ggtitle("Raw \nlog2 counts") + 
      facet_grid(. ~ Treatment  ) +
      ylab("log2(x+1)") + xlab("") +
      theme(title = element_text(size=10), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            plot.title = element_text(hjust = 0.5, face="bold", size=18))

  # SCTtransform
  cat("Reading, summarising and plotting feature counts for:", norm.assay, "assay")
  
  Norm.df <-  GetAssayData(seurat.corrected[[norm.assay]]) %>%
    as.matrix(.) %>% 
    reshape2::melt(., id = rownames(.)) %>% # Convert to long table for easier plotting
    dplyr::rename("id" = 1, "sample" = 2, "value" = 3) %>% 
    tidyr::separate(col ="sample", into = c('Well', 'Treatment', 'Batch'), remove = FALSE) %>% # This take a while due to so many columns created when melting each row to show genes by cell. If only after the average expression. Could probably average before...
    group_by(Batch, Treatment, id) %>% summarise(mean = mean(value))
  ## Update batch names and specify order for plotting
  Norm.df$Batch <- batch.names[Norm.df$Batch]
  Norm.df$Batch <- factor(Norm.df$Batch, level = batch.order)

  ## plot
  Norm.df_plot <- 
    ggplot(data=Norm.df, aes(x=Batch, y=mean, fill=Batch)) + 
      geom_boxplot()  + 
      theme_bw() +
      ggtitle("Normalised \nSCTransformed") + 
      facet_grid(. ~ Treatment  ) +
      ylab("regularized negative binomial regression residuals") + xlab("") +
      theme(title = element_text(size=10), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            plot.title = element_text(hjust = 0.5, face="bold", size=18))


  # Integrated/Corrected
  cat("Reading, summarising and plotting feature counts for:", corrected.assay, "assay")
  
  Corrected.df <-  GetAssayData(seurat.corrected[[corrected.assay]]) %>%
    as.matrix(.) %>% 
    reshape2::melt(., id = rownames(.)) %>% # Convert to long table for easier plotting
    dplyr::rename("id" = 1, "sample" = 2, "value" = 3) %>% 
    tidyr::separate(col ="sample", into = c('Well', 'Treatment', 'Batch'), remove = FALSE) %>% # This take a while due to so many columns created when melting each row to show genes by cell. If only after the average expression. Could probably average before...
    group_by(Batch, Treatment, id) %>% summarise(mean = mean(value))
  ## Update batch names and specify order for plotting
  Corrected.df$Batch <- batch.names[Corrected.df$Batch]
  Corrected.df$Batch <- factor(Corrected.df$Batch, level = batch.order)

  ## plot
  Corrected.df_plot <- 
    ggplot(data=Corrected.df, aes(x=Batch, y=mean, fill=Batch)) + 
      geom_boxplot()  + 
      theme_bw() +
      ggtitle("Seurat batch corrected \nSCTransformed + integrated") + 
      facet_grid(. ~ Treatment  ) +
      ylab("counts \n(batch corrected)") + xlab("") +
      theme(title = element_text(size=10), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face="bold", size=18)) 


  # Combine plots
  cowplot::plot_grid(Raw.df_plot, Norm.df_plot, Corrected.df_plot, ncol =3)
  cat("Saving figure to: ",output.fig)
  ggsave(filename = output.fig, width = 15, height =8)

}
