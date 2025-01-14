##########################################################################
# Purpose: R Script with two functions to heatmap gene expression for 1) correlated genes (csv file provided) or 2) all genes.
# Output: Produces heatmaps and a csv tables that identify which gene belonged to which cluster. The number of heatmaps is optional.
#
# Date: 09.Feb.23
# Version: v.0.1.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## Added option to save images as .svg for HVGs_select_and_PCA()
# To do: 
## 1) Save list of HVGs (genes and TEs)
## 2) Overlay any TEs identified as HVGs on the PCA/ umap plots! (separate function)
## 3) Automate saving of HVGs for k-means by pheatmap and ComplexHeatmap
## 4) Obtain cluster list from ComplexHeatmap like  pheatmap
## 5) Run ClusterProfiler on lists (separate function)

##########################################################################

#Load packages
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(dplyr, data.table, stringr, SingleCellExperiment, scran, pheatmap, ggplot2, cowplot, clustree, tibble)

#################################################################################################################################################
## First function identifies and plots PCA of HVGs and saves list
#################################################################################################################################################

HVGs_select_and_PCA <- function(sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds",
                                chr.to.keep = c("1", "2", "3", "4", "5", "6","7","8", "9","10","11", "12","13","14", "15","16","17", "18","19","20","21", "22", NA), # NA represents TEs 
                                Keep.genes.TEs = c("gene","TE"), # Enables HVG selection to filter for either just gene or TE, or both N.B. Make sure 'chr.to.keep' doesn't filter out TEs if you want to keep them.
                                output.fig.dir = "./figures/colony/HL60/HVGs/",
                                output.file.dir = "./2_results/colony/HL60/HVGs/",
                                output.file.name = 'HVGs_2000.txt',
                                batch.correction.method = "mnnCorrected",
                                PCA.methylation.threshold = 75,
                                Number.of.genes = 2000,
                                size = 3,
                                width = 6,
                                height = 4,
                                alpha = 0.6,
                                svg = FALSE, # If true, images will be saved as both .png and .svg
                                all_cols = c("unt" = "green", "aza" =  "blue", "dac"  = "red",
                                             "rep1" = "orange", "rep2" = "violet", "rep3" = "brown",),
                                set.the.seed = TRUE){

#### Test Variables to delete
#sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds"
#chr.to.keep = c("1", "2", "3", "4", "5", "6","7","8", "9","10","11", "12","13","14", "15","16","17", "18","19","20","21", "22", NA) # NA represents TEs 
#Keep.genes.TEs <- c("gene","TE") # Enables HVG selection to filter for either just gene or TE, or both N.B. Make sure 'chr.to.keep' doesn't filter out TEs if you want to keep them.
#output.fig.dir = "./figures/colony/HL60/HVGs/"
#output.file.dir = "./2_results/colony/HL60/HVGs/"
#output.file.name = 'HVGs_2000.txt'
#batch.correction.method = "mnnCorrected"
#PCA.methylation.threshold = 75
#Number.of.genes = 2000
#size = 3
#width = 6
#height = 4
#alpha = 0.6
#all_cols = c("unt" = "green", "aza" =  "blue", "dac"  = "red",
#             "rep1" = "orange", "rep2" = "violet", "rep3" = "brown")


######################

  # Set seed (for reproducibility of PCA plots)
  if(set.the.seed == TRUE){set.seed(1)}

# Create dir (if required)
  dir.create(output.fig.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output.file.dir, recursive = TRUE, showWarnings = FALSE)

##################################
# Read in file and filter. Rejoin PCA results to main sce.rds file after.
  message("Reading in sce datafile: '", sce.rds.file, "'")
#
  data <- readRDS(sce.rds.file)



## Filter for genes and/or TEs
  cat("\nFiltering for genes and/or TEs: ", Keep.genes.TEs)
  data_filt <- data[rowData(data)$Type %in% Keep.genes.TEs]
## Filter Chromosomes
  chr_to_keep_index <- which(rowData(data_filt)$Chromosome %in% chr.to.keep)
  data_filt <- data_filt[chr_to_keep_index,]
  cat("\nChromosomes kept: \n", unique(rowData(data_filt)$Chromosome))


# Calculate PCA (top X genes). 
  cat("\nCalculating PCA results for top: ", Number.of.genes, "(",Keep.genes.TEs, ")\nUsing the batch corrected dataframe: ", batch.correction.method)
  data_filt_PCA_results <- scater::calculatePCA(data_filt, exprs_values = batch.correction.method, ntop = Number.of.genes)

# Store results back in orignal sce
  New_Dim_Name <- paste0("PCA_",batch.correction.method,"_Top_",Number.of.genes,"_HVGs")
  reducedDims(data)[[New_Dim_Name]] <- data_filt_PCA_results

## Save HVG list from PCA ###############
  PCA_Genes_used <- rownames(attr(data_filt_PCA_results, "rotation")) ## Genes used in PCA
  HVGs_table <- data.frame(HVG_list = PCA_Genes_used)
  HVGs_table$Type <- ifelse(grepl("ENS", HVGs_table$HVG_list), "gene", "TE")


  fwrite(HVGs_table, 
            file = paste0(output.fig.dir, output.file.name), sep = "\t", row.names = F, quote = F)

#########################################
# Plot and save PCA using HVGs

## Add column in sce for above above below methylation threshold
  cat("\nAdding 'high' vs 'low' for methylation level per cell, using the threshold: ", PCA.methylation.threshold, "%\n", sep ="")
  colData(data)$Methylation_level <- ifelse(colData(data)$meanMeth > PCA.methylation.threshold, "High", "Low")

## Modified plotting function ###################################33
    plot_reducedDim <- function(sce, reducedDim = 'PCA', comps = c(1,2), colBy = 'Batch', facet_by = NULL, cols, axis.title = 'PC', ...)
    {
      df <- sce@int_colData$reducedDims[[reducedDim]]
      df <- df[,comps]
      df <- as.data.frame(df)
      colnames(df) <- paste0(axis.title, comps)
      pcs <- colnames(df)
      df$group <- colData(sce)[,colBy]
      df$facet <- colData(sce)[,facet_by]
      cols <- cols[unique(df$group)]
      df$shapes <- colData(sce)[,shapes]
  
      p <- ggplot(df, aes_string(pcs[1], pcs[2])) + geom_point(aes(col = group, shape = shapes),...) +
        scale_colour_manual(values = cols) +
        theme_classic() + labs(col = colBy, shape = shapes.lab)
      if (!is.null(facet_by))
      {
        p <- p + facet_wrap(.~facet, scales = 'free')
      }
      p
    }
###################################################################
## Apply variables to custom PCA function

  output.fig.dir = output.fig.dir
  PCA.name = New_Dim_Name
  shapes = "Methylation_level"
  shapes.lab = paste0("Methylation >",PCA.methylation.threshold,"%")
  size = size
  width = width
  height = height
  alpha = alpha
  all_cols = all_cols
  colBy_options = c('Batch', 'Treatment')
  sce <- data

##
  for (colBy in colBy_options){
    plot_reducedDim(sce = sce, reducedDim = PCA.name,
                    comps = c(1,2), colBy = colBy,
                    cols = all_cols, axis.title = 'PC', size = size, shapes = shapes, alpha = alpha, shapes.lab = shapes.lab)
    filename <- paste0(output.fig.dir, PCA.name, "-By",colBy,".png")

    dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    cat("writing: \n ", filename)
    ggsave(filename =  filename,
           width = width, height = height)
    
    if(svg == TRUE){ # Adds option to save as .svg
      filename <- paste0(output.fig.dir, PCA.name, "-By",colBy,".svg")
      ggsave(filename =  filename,
             width = width, height = height)
    }
  }

  ##############################################################################################################################################
  # Finished
  cat("Finished extracting and plotting PCA for 'top' (",Number.of.genes,") HVGs \nUse the saved HVG table to heatmap, which was stored in: ", paste0(output.fig.dir, output.file.name))
}
################################################################################################################################################





#################################################################################################################################################
# Identify Optimal Kmeans Cluster
#################################################################################################################################################


Kmeans_Optimal_Kmeans_Clusters <- function(sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds",
                                           HVG.table = "./figures/colony/HL60/HVGs/HVGs_2000.txt",
                                           kmeans.elbow = 40, # Number of clusters to test in elbow plots
                                           kmeans.gap.stat = 24, # Number of clusters to test in gap stat plots (takes longer to run than elbow plots)
                                           kmeans.clustree = 9, # Number of clusters to display on clustree. Max currently allowed is 11.
                                           run.Kmean.index = FALSE, # This is set by default to false as it takes a long time.
                                           kmeans.index = 15, # Number of clusters to apply on multiple statistical test to determine 'optimal' number. WARNING - THIS TAKES A LONG TIME!
                                           
                                           output.fig.dir = "./figures/colony/HL60/HVGs/Kmeans_Opt/",
                                           output.file.dir = "./2_results/colony/HL60/HVGs/Kmeans_Opt/",
                                           batch.correction.method = "mnnCorrected",
                                           set.the.seed = TRUE){ # Defaults to true and will set.seed(123) for reproducibility. Set to FALSE to prevent this.){ # Defaults to true and will set.seed(123) for reproducibility. Set to FALSE to prevent this.){}
  
  # TO DELETE TEST VARIABLES ##################
  #sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds"
  #HVG.table = "./figures/colony/HL60/HVGs/HVGs_2000.txt"
  
  #kmeans.elbow = 40 # Number of clusters to test in elbow plots
  #kmeans.gap.stat = 24 # Number of clusters to test in gap stat plots (takes longer to run than elbow plots)
  #kmeans.clustree = 9 # Number of clusters to display on clustree. Max currently allowed is 11.
  #run.Kmean.index = FALSE # This is set by default to false as it takes a long time.
  #kmeans.index = 15 # Number of clusters to apply on multiple statistical test to determine 'optimal' number. WARNING - THIS TAKES A LONG TIME!
  
  #output.fig.dir = "./figures/colony/HL60/HVGs/Kmeans_Opt/TEST/"
  #batch.correction.method = "mnnCorrected"
  #set.the.seed = TRUE # Defaults to true and will set.seed(123) for reproducibility. Set to FALSE to prevent this.
  
  ############################################
  
  # set.seed if desired
  if(set.the.seed == TRUE){set.seed(123)}
  
  # Create dir (if required)
  dir.create(output.fig.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output.file.dir, recursive = TRUE, showWarnings = FALSE)
  
  # Set seed (for reproducibility of PCA plots)
  set.seed(1)
  
  # Import HVG table 
  HVG_table <- fread(HVG.table)
  HVG_list <- HVG_table %>% .[[1]]
  
  # Import sce
  message("Reading in sce datafile: '", sce.rds.file, "'")
  data <- readRDS(sce.rds.file)
  
  ## Filter for batch correction and HVG list
  message("\nUsing the assay/batch corrected data: ",batch.correction.method)
  Gene_Matrix <- as.data.frame(assay(data, batch.correction.method)) # Filter for batch correction method
  Gene_Matrix <- Gene_Matrix[row.names(Gene_Matrix) %in% HVG_list,] # Filter for genes
  message("\nAfter filtering, retained and performing heatmapping on ",nrow(Gene_Matrix), " HVGs")
  
  # Prepare matrices
  ## Apply row normalisation - it's already per batch normalised and batch corrected. So need to perform vst normalisation.
  message("Applying row normalisation by: \n1) Mean centering \n2) Z-scaling")
  ###  1) Row normalised by Mean centering (x-RowMean)
  Mean_centred <- Gene_Matrix - rowMeans(Gene_Matrix) 
  ###  2) Row normalised by Z score (x-Rowmean/ SD)
  Z_scored <- t(scale(t(Gene_Matrix)))
  
  cat("\nProviding quick summary of optimal Kmeans\n")
  library(factoextra)
  
  ### ELBOW PLOTS ###################################################################################
  message("\nCalculating elbow plots \nMax kmeans being tested: ", kmeans.elbow)
  
  elbow_raw <- fviz_nbclust(t(Gene_Matrix), kmeans, method = "wss", k.max = kmeans.elbow) + theme_minimal() + ggtitle("'Raw' counts")
  elbow_mean <- fviz_nbclust(t(Mean_centred), kmeans, method = "wss", k.max = kmeans.elbow) + theme_minimal() + ggtitle("Mean centred")
  elbow_zscore <- fviz_nbclust(t(Z_scored), kmeans, method = "wss", k.max = kmeans.elbow) + theme_minimal() + ggtitle("z-score")
  
  
  elbow_plot_row <- plot_grid(elbow_raw,
                              elbow_mean,
                              elbow_zscore, label_size = 12, nrow = 1) 
  
  title <- ggdraw() + draw_label(
    "Elbow plots of Kmeans clustering",
    fontface = 'bold',
    x = 0, hjust = 0,size = 24) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))
  
  plot_grid(title,
            elbow_plot_row, 
            label_size = 12, ncol = 1,rel_heights = c(0.1, 1)) + 
    theme(plot.background = element_rect(fill= "white", color = NA))
  
  
  ggsave(filename =  paste0(output.fig.dir,"Kmeans_Elbow_plots.png"), height = 8, width = 24)
  
  
  ### GAP STATISTIC ###################################################################################
  message("\nCalculating Gap statistic\nMax kmeans being tested: ", kmeans.gap.stat)
  message("\nThis can take a little while... If it's taking too long. Reduce the number kmeans clusters (or go into the script and change number of Bootstrapps 'B'.\n")
  
  library(cluster)
  gap_stat_raw <- cluster::clusGap(t(Gene_Matrix), FUN = kmeans, nstart = 20, K.max = kmeans.gap.stat, B = 30)
  gap_stat_mean <- cluster::clusGap(t(Mean_centred), FUN = kmeans, nstart = 20, K.max = kmeans.gap.stat, B = 30)
  gap_stat_zscore <- cluster::clusGap(t(Z_scored), FUN = kmeans, nstart = 20, K.max = kmeans.gap.stat, B = 30)
  
  gap_stat_raw_plot <- fviz_gap_stat(gap_stat_raw) + theme_minimal() + ggtitle("'Raw' counts")
  gap_stat_mean_plot <- fviz_gap_stat(gap_stat_mean) + theme_minimal() + ggtitle("Mean centred")
  gap_stat_zscore_plot <- fviz_gap_stat(gap_stat_zscore) + theme_minimal() + ggtitle("Z-score")
  
  gap_plot_row <- plot_grid(gap_stat_raw_plot,
                            gap_stat_mean_plot,
                            gap_stat_zscore_plot, label_size = 12, nrow = 1) 
  
  title <- ggdraw() + draw_label(
    "fviz_gap_stat: Gap Statistic plots of Kmeans clustering",
    fontface = 'bold',
    x = 0, hjust = 0,size = 24) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7))
  
  plot_grid(title,
            gap_plot_row, 
            label_size = 12, ncol = 1,rel_heights = c(0.1, 1)) + 
    theme(plot.background = element_rect(fill= "white", color = NA))
  
  ggsave(filename =  paste0(output.fig.dir,"Kmeans_Gap_Stat__plots.png"), height = 8, width = 24)
  
  
  ### CLUSTREE ###################################################################################
  message("\nCalculating clustree Kmeans PCA \nMax kmeans being plotted: ", kmeans.clustree)
  
  for(matrix in 1:2 ){
    if(matrix ==1 ){
      to_test_clustree <- as.data.frame(t(Mean_centred))
    } else if(matrix ==2){
      to_test_clustree <- as.data.frame(t(Z_scored))
    }
    
    
    
    tmp <- NULL
    k_max_to_test <- 11
    
    for (k in 1:k_max_to_test){
      tmp[k] <- kmeans(to_test_clustree, k, nstart = 30)
    }
    
    df <- data.frame(tmp)
    # add a prefix to the column names
    colnames(df) <- seq(1:k_max_to_test)
    colnames(df) <- paste0("k",colnames(df))
    # get individual PCA
    df.pca <- prcomp(df, center = TRUE, scale. = FALSE) 
    # df.pca <- prcomp(t(to_test_clustree), center = TRUE, scale. = FALSE)
    ind.coord <- df.pca$x
    # ind.coord <- df.pca$rotation
    # ind.coord <- reducedDim(data2, "PCA")[,1:2]
    ind.coord <- ind.coord[,1:2]
    df <- bind_cols(as.data.frame(df), as.data.frame(ind.coord))
    clustree(df, prefix = "k")
    
    kmean_subset <- 1:kmeans.clustree
    
    df_subset <- df %>% dplyr::select(kmean_subset,(ncol(df)-1):ncol(df))
    PCA_clustree <- clustree_overlay(df_subset, prefix = "k", x_value = "PC1", y_value = "PC2")
    
    library(cowplot)
    overlay_list <- clustree_overlay(df_subset, prefix = "k", x_value = "PC1",
                                     y_value = "PC2", plot_sides = TRUE)
    
    PC1_clustree <- overlay_list$x_side
    
    PC2_clustree <- overlay_list$y_side
    
    bottom_row <- plot_grid(PC1_clustree + guides(fill = "none", col = "none"), 
                            PC2_clustree  + guides(fill = "none", col = "none"), 
                            labels = c('PC1', 'PC2'), label_size = 12)
    
    
    
    plot_grid(PCA_clustree + guides(fill = "none"), 
              bottom_row, 
              labels = c('PCA clustree', ''), label_size = 12, ncol = 1) + 
      theme(plot.background = element_rect(fill= "white", color = NA))
    
    if(matrix ==1 ){
      ggsave(filename =  paste0(output.fig.dir,"Kmeans_Clustree_assess_by_PCA_Mean.png"), height = 14, width = 10)
    } else if(matrix ==2){
      ggsave(filename =  paste0(output.fig.dir,"Kmeans_Clustree_assess_by_PCA_Zscore.png"), height = 14, width = 10)
    }
    
    
  }
  
  ### Combined index score for optimal Kmeans ###################################################################################
  
  if(run.Kmean.index == TRUE){
    cat("Calculating and combining kmeans index score across multiple test as 'run.Kmean.index == TRUE'")
    cat("\nThe number of kmean clusters being assessed is automatically set to compare 2 up to 15. And currently only assess the scaled (z-score) data.")
    message("\nThis can take some time to run... The other plots should already be completed and saved for your inspection while this finishes.\n")
    # to_test_clustree <- as.data.frame(t(Mean_centred))
    to_test_clustree <- as.data.frame(t(Z_scored))
    
    selected <- c( "kl", "ch", "hartigan",  "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "frey", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw") 
    results <- NULL
    for (i in 1:length(selected)) {
      results[[i]] <- try(NbClust(to_test_clustree, min.nc=2, max.nc=15, method="kmeans", index=selected[i]))
    }
    
    score.index <- NULL
    for(i in 1:length(results) ){
      score.index <- c(score.index, results[[i]]$Best.nc[[1]])
    }
    
    png(paste0(output.fig.dir,"Kmeans_optimal_Kmeans_cluster_NbClust_Zscore.png"), width = 800, height = 600)
    hist(score.index,
         main="Frequency of 'optimal' Kmeans clusters from NbClust indexes",
         xlab="Kmeans cluster number",
         col="darkmagenta", breaks = 15, xlim = c(1,15),
         freq=TRUE) 
    dev.off()
  }
  ###############################################################################################################################################  
  cat("\nFinished all plots. Now check the figures to choose an 'optimal' number of Kmean cluster to apply.")
  cat("\nFigures are located in:", output.fig.dir,"\n")
  
}
#################################################################################################################################################






#################################################################################################################################################
# Heatmap using gene list
#################################################################################################################################################


Heatmap_from_list <- function(sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds",
                              HVG.table = "./figures/colony/HL60/HVGs/HVGs_2000.txt",
                              methylation.threshold = 75,
                              kmeans.clusters = 8,
                              kmeans.iter = 50,
                              kmeans.nstart = 30,
                              output.fig.dir = "./figures/colony/HL60/HVGs/",
                              output.file.dir = "./2_results/colony/HL60/HVGs/",
                              batch.correction.method = "mnnCorrected",
                              plot.Heatmap.Row.Normalisations = c("Zscored", "Mean_Centred", "Raw"), # Can be restricted to just one or two options. c("Zscored", "Mean_Centred", "Raw")
                              plot.Pheatmap.style = TRUE,
                              plot.Complex.Heatmap.NoRowCluster.style = TRUE,
                              plot.Complex.Heatmap.WithRowCluster.style = TRUE,
                              set.the.seed = TRUE){ # Defaults to true and will set.seed(123) for reproducibility. Set to FALSE to prevent this.){}

# TO DELETE TEST VARIABLES ##################
#sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds"
#HVG.table = "./figures/colony/HL60/HVGs/HVGs_2000.txt"

#methylation.threshold = 75
#kmeans.clusters = 8
#kmeans.iter = 50
#kmeans.nstart = 30

#output.fig.dir = "./figures/colony/HL60/HVGs/"
#output.file.dir = "./2_results/colony/HL60/HVGs/"
#batch.correction.method = "mnnCorrected"
#set.the.seed = TRUE # Defaults to true and will set.seed(123) for reproducibility. Set to FALSE to prevent this.

############################################

# set.seed if desired (for reproducibility of PCA plots)
  if(set.the.seed == TRUE){set.seed(123)}

############################################

# Create dir (if required)
  dir.create(output.fig.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output.file.dir, recursive = TRUE, showWarnings = FALSE)


# Import HVG table 
  HVG_table <- fread(HVG.table)
  HVG_list <- HVG_table %>% .[[1]]

# Import sce
  message("Reading in sce datafile: '", sce.rds.file, "'")
  data <- readRDS(sce.rds.file)

## Filter for batch correction and HVG list
  message("\nUsing the assay/batch corrected data: ",batch.correction.method)
  data_filt <- data[rownames(data) %in% HVG_list] # Filter for genes
  Gene_Matrix <- as.data.frame(assay(data_filt, batch.correction.method)) # Filter for batch correction method
  message("\nAfter filtering, retained and performing heatmapping on ",nrow(Gene_Matrix), " HVGs")


  
  
# Prepare matrices

## Using colData as samples sheet, for most up to date QC.
### Create high vs low methlyation level column (defined methylation level). 
  cat("\nAdding 'high' vs 'low' for methylation level per cell, using the threshold: ", methylation.threshold, "%\n", sep ="")
  colData(data_filt)$Methylation_level <- ifelse(colData(data_filt)$meanMeth > methylation.threshold, "High", "Low")

  colData_filtered <- as.data.frame(colData(data_filt)) %>%
    arrange(desc(Treatment), desc(meanMeth)) %>%
    dplyr::select(Treatment, Methylation_level, meanMeth, batch) 

  message("\nTotal sample numbers (post-QC) in each treatment and methylation level:")
  print(colData_filtered %>% 
        group_by(Treatment, Methylation_level) %>% 
        tally())

## Re-order matrix to be in order of treatment and descending methylation level!
  Gene_Matrix <- as.data.frame(Gene_Matrix) %>% dplyr::select(row.names(colData_filtered))

## Apply row normalisation - it's already per batch normalised and batch corrected. So need to perform vst normalisation.
  message("Applying row normalisation by: \n1) Mean centering \n2) Z-scaling")
###  1) Row normalised by Mean centering (x-RowMean)
  Mean_centred <- Gene_Matrix - rowMeans(Gene_Matrix) 
###  2) Row normalised by Z score (x-Rowmean/ SD)
  Z_scored <- t(scale(t(Gene_Matrix)))


  ### SETTING UP WITH Z_scored only first. Will loop for other ones after
  
  for(HVG_matrix in plot.Heatmap.Row.Normalisations){ # Looping 'raw', 'mean centred' and 'z scored' count tables for clustering and heatmapping. Can restrict to one or all three - c("Zscored", "Mean_Centred", "Raw")
    #HVG_matrix_to_run <- HVG_matrix
    if(HVG_matrix == "Raw"){HVG_matrix_to_run <- Gene_Matrix
      } else if(HVG_matrix == "Mean_Centred"){HVG_matrix_to_run <- Mean_centred
      } else if(HVG_matrix == "Zscored"){HVG_matrix_to_run <- Z_scored
    }

  
    # Create dir for heatmap being plotted
    output.fig.dir.heatmap <-paste0(output.fig.dir,HVG_matrix,"/")
    dir.create(output.fig.dir.heatmap, recursive = TRUE, showWarnings = FALSE)
  

# Perform K-means clustering
    cat("Performing Kmeans cluster with:\ncluster = ", kmeans.clusters, "\niterations = ",kmeans.iter, "\nnstart =", kmeans.nstart)
    kmeans_run <- kmeans(HVG_matrix_to_run, centers = kmeans.clusters, iter.max = kmeans.iter, nstart = kmeans.nstart)$cluster
  
   kmeans_run_df <-  as.data.frame(kmeans_run)  %>% 
      mutate(gene = row.names(.)) %>% 
      dplyr::select(gene, cluster = kmeans_run) %>% 
      arrange(cluster)
  
## Save kmeans cluster table
    cat("\nSaving kmeans cluster list before running heatmaps. Automatically saved: \n",paste0(output.fig.dir.heatmap,"Kmeans_cluster_list_",HVG_matrix,".txt"))
    fwrite(kmeans_run_df, file = paste0(output.fig.dir.heatmap,"Kmeans_cluster_list_",HVG_matrix,".txt"), sep = "\t", row.names = F)
  
## Reorder matrix rows by kmeans clusters
    row_list_to_order_by <- kmeans_run_df$gene
    HVG_matrix_to_run_kmeans_ordered <- HVG_matrix_to_run[row_list_to_order_by, ]
  


# Prepare heatmap info (annotations and colours)

## Row anno 
  # Need to re-order sce row order to match kmeans order and add cluster column
    annotation_row <- rowData(data_filt)[row_list_to_order_by,] %>% 
      as.data.frame(.) %>% 
      dplyr::select(Type)
    annotation_row <- merge(as.data.frame(kmeans_run_df), annotation_row, by=0, all=TRUE, sort =FALSE) %>% ## Adding cluster label
      column_to_rownames(., var = "Row.names") %>% 
      dplyr::select(-gene)
    annotation_row$cluster <- factor(annotation_row$cluster, levels=1:kmeans.clusters)
  # annotation_row <- data.frame(Type = rowData(data_filt)$Type)
  # row.names(annotation_row) <- row.names(data_filt)

## Column anno - Already sorted above when arranging row order
    colData_filtered
    annotation_column <- colData_filtered

## Colours of annos
  ### Annotation colours
    library(RColorBrewer)
  
    Cluster_colours <- brewer.pal(n = kmeans.clusters, name = "Dark2") ## For the colour of the clusters 
    names(Cluster_colours) <- as.character(1:kmeans.clusters)
    anno_colors = list(
      Treatment= c(dac = "red", unt = "green", aza = "blue"),
      Methylation_level = c(Low = "black", High = "orchid"),
      #Type = c(gene = "black"),
      batch =  c("rep1" = "orange", "rep2" = "violet", "rep3" = "brown"),
      Type = c(gene = "azure4", TE = "blue"),
      meanMeth = c("light blue", "darkgreen"),
      cluster = Cluster_colours)
  
  ### Matrix Colours
    Matrix_to_colour_scale <- HVG_matrix_to_run_kmeans_ordered
    paletteLength <- 200
    myColor <- colorRampPalette(c("blue", "light blue", "white", "orange", "red"))(paletteLength)
    myBreaks <- c(seq(min(Matrix_to_colour_scale), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(Matrix_to_colour_scale)/paletteLength, max(Matrix_to_colour_scale), length.out=floor(paletteLength/2)))



# Plot heatmaps
  
  
  ## PHEATMAP
    if(plot.Pheatmap.style == TRUE){
      PH_Kmeans_no_RowCluster <- pheatmap::pheatmap(HVG_matrix_to_run_kmeans_ordered,
                     # kmeans_k = 8, # decide how many clusters (set to NA if not concerned)
                     #treeheight_row = 1,
                     color = myColor,
                     breaks = myBreaks,
                     
                     main = paste0("Heatmap of top 2000 HVGs (", HVG_matrix,")"),
                     show_rownames = F, show_colnames = F,
                     cutree_cols = 1,
                     clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",  scale = "none", cluster_cols = F, cluster_rows = F,
                     annotation_row = annotation_row,
                     annotation_col = annotation_column,
                     annotation_colors = anno_colors,
                     #color = colorRampPalette(c("blue", "light blue", "white", "light orange", "orange", "red"))(30),
                     
                     annotation_names_row = F,
                     fontsize = 8.5)
   
    ### Save
      ggsave(plot=PH_Kmeans_no_RowCluster,
            file=paste0(output.fig.dir.heatmap,"P.Heatmap_HVGs_Kmeans_",HVG_matrix,".svg"), 
            width=8, height=9)
    }
  ## Complex Heatmap - pheatmap style
    if(plot.Complex.Heatmap.NoRowCluster.style == TRUE){
      CH_Kmeans_no_RowCluster <- ComplexHeatmap::pheatmap(as.matrix(HVG_matrix_to_run_kmeans_ordered),
                           #row_km = 10, # decide how many clusters (set to NA if not concerned)
                           #row_km_repeats = 200, # Will repeat this many times and may reduce number of k-mean clusters depending on result...
                           #treeheight_row = 1,
                           color = myColor,
                           breaks = myBreaks,
                           
                           main = paste0("Heatmap of top 2000 HVGs (", HVG_matrix,")"),
                           show_rownames = F, show_colnames = F,
                           #cutree_cols = 1,
                           #clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",  
                           scale = "none", cluster_cols = F, cluster_rows = F,
                           annotation_row = annotation_row,
                           annotation_col = annotation_column,
                           annotation_colors = anno_colors,
                           #color = colorRampPalette(c("blue", "light blue", "white", "light orange", "orange", "red"))(30),
                           
                           annotation_names_row = F,
                           fontsize = 8.5)
      ### Save
      png(filename= paste0(output.fig.dir.heatmap,"Complex.Heatmap_HVGs_Kmeans_",HVG_matrix,".png"), width=9,height=10,units="in",res=1200)
        ComplexHeatmap::draw(CH_Kmeans_no_RowCluster)
      dev.off()
      svg(filename= paste0(output.fig.dir.heatmap,"Complex.Heatmap_HVGs_Kmeans_",HVG_matrix,".svg"), width=9,height=10)
        ComplexHeatmap::draw(CH_Kmeans_no_RowCluster)
      dev.off()
      pdf(file= paste0(output.fig.dir.heatmap,"Complex.Heatmap_HVGs_Kmeans_",HVG_matrix,".pdf"), width=9,height=10)
        ComplexHeatmap::draw(CH_Kmeans_no_RowCluster)
      dev.off()
    }


    ## Complex Heatmap - Clustering within Kmeans clusters
    if(plot.Complex.Heatmap.WithRowCluster.style == TRUE){
      CH_Kmeans_cluster_within <- ComplexHeatmap::pheatmap(as.matrix(HVG_matrix_to_run_kmeans_ordered),
                           #row_km = 10, # decide how many clusters (set to NA if not concerned)
                           #row_km_repeats = 200, # Will repeat this many times and may reduce number of k-mean clusters depending on result...
                           #treeheight_row = 1,
                           color = myColor,
                           breaks = myBreaks,
                           
                           main = paste0("Heatmap of top 2000 HVGs (", HVG_matrix,") with clustering within kmeans clusters"),
                           show_rownames = F, show_colnames = F,
                           #cutree_cols = 1,
                           #clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",  
                           scale = "none", cluster_cols = F, 
                           row_split = annotation_row$cluster, cluster_row_slices = T, 
                           #cluster_rows    = F, # Should stop it reordering the annotation groups....
                           show_row_dend = F,
                           
                           annotation_row = annotation_row,
                           annotation_col = annotation_column,
                           annotation_colors = anno_colors,
                           
                           annotation_names_row = F,
                           fontsize = 8.5)

        ### Save
        png(filename= paste0(output.fig.dir.heatmap,"Complex.Heatmap_HVGs_Kmeans_within_clusters_",HVG_matrix,".png"), width=9,height=10,units="in",res=1200)
          ComplexHeatmap::draw(CH_Kmeans_cluster_within)
        dev.off()
        svg(filename= paste0(output.fig.dir.heatmap,"Complex.Heatmap_HVGs_Kmeans_within_clusters_",HVG_matrix,".svg"), width=9,height=10)
          ComplexHeatmap::draw(CH_Kmeans_cluster_within)
        dev.off()
        pdf(file= paste0(output.fig.dir.heatmap,"Complex.Heatmap_HVGs_Kmeans_within_clusters_",HVG_matrix,".pdf"), width=9,height=10)
          ComplexHeatmap::draw(CH_Kmeans_cluster_within)
        dev.off()
    }   
        ## Save params used
        sink(paste0(output.fig.dir.heatmap, "Kmeans_cluster_list_",HVG_matrix,"_params_used.txt"))
          L1 <- cat("Params used... (", HVG_matrix,")", sep = "")
          L2 <- cat("\nHVG table:", HVG.table, sep = "\t")
          L3 <- cat("\nbatch.correction.method:", batch.correction.method, sep = "\t")
          L4 <- cat("\nNumber of Kmeans clusters:", kmeans.clusters, sep = "\t")
          L5 <- cat("\nkmeans max iterations:", kmeans.iter, sep = "\t")
          L6 <- cat("\nKmeans.nstart:", kmeans.nstart, sep = "\t")
        sink()
    }
}
  
#############################################################################################################################################################################################