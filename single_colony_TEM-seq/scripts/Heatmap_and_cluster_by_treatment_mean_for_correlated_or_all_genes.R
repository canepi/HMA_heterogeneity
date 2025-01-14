##########################################################################
# Purpose: R Script with two functions to heatmap gene expression averaged across treatment groups, for 1) correlated genes (csv file provided) or 2) all genes.
# Output: Produces heatmaps and a csv tables that identify which gene belonged to which cluster. The number of heatmaps is optional.
#
# Date: 23.Jan.23
# Version: v.0.0.2
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) This relies on the sce.rds already being QC filtered. 
## 2) Cleaned and added set.seed option
# To do: 
##########################################################################

## Packages 
library(dplyr)
#library(plotly)
library(data.table)
library(stringr)
library(SingleCellExperiment)
library(scran)
library(pheatmap)

###########################################################################################################################################################################################################


Heatmap_by_treatment_group_corr_genes <- function(sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds",
                                                  correlated.gene.results = "./2_results/colony/HL60/DAC_high_vs_low/Methylation_Exp_mnnCorrected_correlation.csv",
                                                  padj.threshold = 0.05,
                                                  corr.est.threshold = 0.4,
                                                  methylation.threshold = 70,
                                                  kmeansk.chosen = 4,
                                                  # sample.sheet = "1_data/colony/HL60/Sample_analysis_sheet.csv",
                                                  treatments = "dac",
                                                  output.fig.dir = "./figures/colony/HL60/DAC_high_vs_low/",
                                                  output.file.dir = "./2_results/colony/HL60/DAC_high_vs_low/",
                                                  batch.correction.method = "mnnCorrected",
                                                  Number.of.plots = 5,
                                                  set.the.seed = TRUE){ # Each run produces slight variations in the grouping of genes, but starts to settle/repeat ~5 runs. This may differ depending on your dataset size and potentially signal noise.
# Variables to delete

#sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds"
#correlated.gene.results = "./2_results/colony/HL60/DAC_high_vs_low/Methylation_Exp_mnnCorrected_correlation.csv"
#padj.threshold = 0.05
#corr.est.threshold = 0.4
#methylation.threshold = 70
#sample.sheet = "1_data/colony/HL60/Sample_analysis_sheet.csv"
#treatments = "dac"
#output.fig.dir = "./figures/colony/HL60/DAC_high_vs_low/"
#output.file.dir = "./2_results/colony/HL60/DAC_high_vs_low/"
#batch.correction.method = "mnnCorrected"


  # Set seed (for reproducibility of PCA plots)
  if(set.the.seed == TRUE){set.seed(123)}

  ###### Create figure folder ########

  mkdirs <- function(fp) {
    if(!file.exists(fp)) {
      mkdirs(dirname(fp))
      dir.create(fp)
    }
  } 
  mkdirs(output.fig.dir)
  mkdirs(output.file.dir)

  ##################################
  # Read in correlated results, then filter for significant hits and split into correlation direction (and combined)
  message("Using the correlated results table: '", correlated.gene.results, "'")
  message("\nFiltering thresholds chosen:\nPadj <= ", padj.threshold, "\nand", "\ncorr.est >= ", corr.est.threshold, " OR corr.est <= -", corr.est.threshold)
  #
  y.updated <- read.csv(correlated.gene.results)

  Corr.value.tested <- corr.est.threshold
  padj.threshold <- padj.threshold
  Pos.corr.cells <- y.updated %>% filter(p.adj <= padj.threshold & cor.val.estimate >= Corr.value.tested) %>% mutate(Gene_corr_group = "Pos_corr_cells")
  neg.corr.cells <- y.updated %>% filter(p.adj <= padj.threshold & cor.val.estimate <= -Corr.value.tested) %>% mutate(Gene_corr_group = "Neg_corr_dac_cells")
  both.corr.cells <- y.updated %>% filter(p.adj <= padj.threshold & cor.val.estimate >= Corr.value.tested | cor.val.estimate <= -Corr.value.tested) %>% mutate(Gene_corr_group = "Pos.Neg_corr_cells")

  ##########################################################
  # Import gene list, filtering for genes positively and negatively correlated AND batch correction method  
  message("Reading in sce datafile: '", sce.rds.file, "'")
  #
  data <- readRDS(sce.rds.file)
  # Filter genes for chosen set identified (by correlation analysis)
  message("Filtering for genes/TEs meeting the chosen correlation thresholds, which totalled: ", length(both.corr.cells$ensembl.TE))
  #
  To_filter <- both.corr.cells # Tested on Pos.corr.dac.cells
  data_filt <- data[rownames(data) %in% To_filter$ensembl.TE, ]

  # Filter for batch correction method
  message("\nand using the assay/batch corrected data: ",batch.correction.method)
  #
  batch_correction_method <- batch.correction.method
  Gene_Matrix <-   assay(data_filt, batch_correction_method)


  # Using colData as samples sheet, for most up to date QC.
  # Create high vs low methlyation level column (defined methylation level). 
  message("Defining 'high' methylated samples as methylation > ", methylation.threshold, "%")
  #
  colData_filtered <- as.data.frame(colData(data)) %>% 
    mutate(Methylation_level = ifelse(meanMeth > methylation.threshold, "high","low")) %>%
    arrange(desc(Treatment), desc(Methylation_level)) %>%
    dplyr::select(Treatment, Methylation_level) 


  message("Total sample numbers (post-QC) in each treatment and methylation level:")
  print(
    colData_filtered %>% group_by(Treatment, Methylation_level) %>% tally()
    )

  # Should save this table

  #################################################################################################
  # Filter sample IDs for high vs low methylation levels in treatment groups
  message("Filtering into treatment groups, and high vs. low methylation dac.")
  #
  unt <-  colData_filtered %>% filter(Treatment == "unt") %>% row.names(.)
  dac_high  <- colData_filtered %>% filter(Treatment == "dac" & Methylation_level == "high") %>% row.names(.)
  dac_low  <- colData_filtered %>% filter(Treatment == "dac" & Methylation_level == "low") %>% row.names(.)
  aza  <- colData_filtered %>% filter(Treatment == "aza") %>% row.names(.)
  #aza_high  <- colData_filtered %>% filter(Treatment_group == "aza" & Methylation_level == "high") %>% row.names(.)
  #aza_low  <- colData_filtered %>% filter(Treatment_group == "aza" & Methylation_level == "low") %>% row.names(.)


  # Matrix of average values for each treatment group (as identified above)
  message("Creating count matix to heatmap, based oon average count for each treatment group (and high vs low dac).")
  #
  Mean_Matrix <- data.frame(unt =rowMeans(Gene_Matrix[,unt], na.rm = TRUE),
                            dac_high =rowMeans(Gene_Matrix[,dac_high], na.rm = TRUE),
                            dac_low =rowMeans(Gene_Matrix[,dac_low], na.rm = TRUE),
                            aza =rowMeans(Gene_Matrix[,aza], na.rm = TRUE)
  )

  # 2) Row normalised by Mean centering (x -RowMean)
  message("Mean centering data")
  #
  Mean_Matrix_Mean_centred <- Mean_Matrix - rowMeans(Mean_Matrix) 


  message("Preparing heatmapping info...")

  ## Group Info
  Mean_matrix_sample_sheet <- data.frame(Treatment = c("unt", "dac","dac","aza"),
                                       Methylation_level = c("high", "high", "low", "low"))
  row.names(Mean_matrix_sample_sheet) <- colnames(Mean_Matrix)

  # Row annotations
  annotation_row <- data.frame(Type = rowData(data_filt)$Type)
  row.names(annotation_row) <- row.names(data_filt) 


  # Set colours
  Matrix_to_colour_scale <- Mean_Matrix
  paletteLength <- 200
  myColor <- colorRampPalette(c("blue", "light blue", "white", "orange", "red"))(paletteLength)
  myBreaks <- c(seq(min(Matrix_to_colour_scale), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(Matrix_to_colour_scale)/paletteLength, max(Matrix_to_colour_scale), length.out=floor(paletteLength/2)))

  Matrix_to_colour_scale <- Mean_Matrix_Mean_centred
  myBreaks_mean_centred <- c(seq(min(Matrix_to_colour_scale), 0, length.out=ceiling(paletteLength/2) + 1), 
                           seq(max(Matrix_to_colour_scale)/paletteLength, max(Matrix_to_colour_scale), length.out=floor(paletteLength/2)))

  ann_colors = list(
    Treatment_group = c(dac = "red", unt = "green", aza = "blue"),
    Methylation_level = c(low = "black", high = "orchid"),
    Type = c(gene = "green", TE = "blue"))



  ### Run/loop this five times, so the user can see how the clustering algorithm may slightly alter groups upon each run. Seems to repeat every 4-5 runs.
  message("Plotting (and saving results) from 5 repeated heatmaps as clustering lustering algorithm may slightly alter groups upon each run. \nThis seems to stabalised every 4-5 runs.")
  message("\'Stabalisation'/variation seems to depend on the number of genes being plotted AND number of clusters to be defined.")
  message("\nKmeans (number of groups/clusters) has been set as: ", kmeansk.chosen)

  for(Run in 1:5){

    Mean_Matrix_rowMeans_with_kmeans <- 
      pheatmap(Mean_Matrix_Mean_centred,
               kmeans_k = kmeansk.chosen, # decide how many clusters (set to NA if not concerned)
               #treeheight_row = 1,
               color = myColor,
               breaks = myBreaks_mean_centred,
             
               main = paste0(batch.correction.method," - genes Pos. and Neg. Corr. with Methylation in DAC cells \n showing average expression in cell treatment groups (Mean Centred)"),
               show_rownames = T, show_colnames = F,
               cutree_cols = 1,
               clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",  scale = "none", cluster_cols = F,
               annotation_row = annotation_row,
               annotation_col = Mean_matrix_sample_sheet,
               annotation_colors = ann_colors,
               #color = colorRampPalette(c("blue", "light blue", "white", "light orange", "orange", "red"))(30),
             
               annotation_names_row = F,
               fontsize = 8.5)
               #filename = paste0(output.fig.dir,"Heatmap_DAC_Pos.Neg_",Corr.value.tested,"_",batch_correction_method,"_Average_across_treatment_groups_MeanCentred_4Km_Run", Run,".png"))
  
  
    ggsave(plot=Mean_Matrix_rowMeans_with_kmeans,
           file=paste0(output.fig.dir,"Heatmap_DAC_Pos.Neg_",Corr.value.tested,"_",batch_correction_method,"_Average_across_treatment_groups_MeanCentred_4Km_Run", Run,".svg"), 
           width=8, height=9)
    ggsave(plot=Mean_Matrix_rowMeans_with_kmeans,
           file=paste0(output.fig.dir,"Heatmap_DAC_Pos.Neg_",Corr.value.tested,"_",batch_correction_method,"_Average_across_treatment_groups_MeanCentred_4Km_Run", Run,".png"), 
           width=8, height=9)

  ### Extract clusters from heatmap ###########################################################################
  # Extract cluster number for each genes and update original correlation table
  ## This bit is more specific to the NHMRC grant

  # Set which pheatmap to extract cluster from
    Pheatmap_to_explore <- Mean_Matrix_rowMeans_with_kmeans
  # Obtain Matrix of clusters
    Gene_Clusters <- as.data.frame(Pheatmap_to_explore$kmeans$cluster) %>% 
      tibble::rownames_to_column( "genes") %>%
      dplyr::rename("cluster_number" = "Pheatmap_to_explore$kmeans$cluster")

#colnames(Gene_Clusters) <- c("genes","cluster_number")
#write.csv(Gene_Clusters,
#          paste0(output.file.dir,"Gene_Clusters_from_Heatmap-Pos.Neg_",Corr.value.tested,"_",batch_correction_method,"_Average_across_treatment_groups_MeanCentred_4.csv"),
 #         row.names = F)

  # Update correlation  results with gene clusters
    y.updated2 <- 
      y.updated %>%
      dplyr::left_join(Gene_Clusters, by = c("ensembl.TE" = "genes"))

# Save updated file
    write.csv(y.updated2, 
              file = paste0(fs::path_ext_remove(correlated.gene.results), "_and_heatmap_cluster_Run",Run,".csv"),
              row.names = F)

  message("Finished plotting heatmap and saved the gene clusters for run ", Run, "\nby updating correlation table and storing as:\n'", fs::path_ext_remove(correlated.gene.results), "_and_heatmap_cluster_Run",Run,".csv", "'")
  }
message("Finished plotting five different heatmaps, check each one to see clustering consistency (of genes) and which heatmap + cluster to investigate further")
}





###########################################################################################################################################################################################################

Heatmap_by_treatment_group_all_genes <- function(sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds",
                                                 #correlated.gene.results = "./2_results/colony/HL60/DAC_high_vs_low/Methylation_Exp_mnnCorrected_correlation.csv",
                                                 #padj.threshold = 0.05,
                                                 #corr.est.threshold = 0.4,
                                                 methylation.threshold = 70,
                                                 kmeansk.chosen = 8,
                                                 # sample.sheet = "1_data/colony/HL60/Sample_analysis_sheet.csv",
                                                 treatments = "dac",
                                                 output.fig.dir = "./figures/colony/HL60/DAC_high_vs_low/all_genes/",
                                                 output.file.dir = "./2_results/colony/HL60/DAC_high_vs_low/all_genes/",
                                                 batch.correction.method = "mnnCorrected",
                                                 Number.of.plots = 5,
                                                 set.the.seed = TRUE){
  
  




# To delete tester variables
#sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds"
#correlated.gene.results = "./2_results/colony/HL60/DAC_high_vs_low/Methylation_Exp_mnnCorrected_correlation.csv"
#padj.threshold = 0.05
#corr.est.threshold = 0.4
#methylation.threshold = 70
#kmeansk.chosen = 8
# sample.sheet = "1_data/colony/HL60/Sample_analysis_sheet.csv",
#treatments = "dac"
#output.fig.dir = "./figures/colony/HL60/DAC_high_vs_low/all_genes/"
#output.file.dir = "./2_results/colony/HL60/DAC_high_vs_low/all_genes/"
#batch.correction.method = "mnnCorrected"

  # Set seed (for reproducibility of PCA plots)
  if(set.the.seed == TRUE){set.seed(123)}

  ###### Create figure folder ########

  mkdirs <- function(fp) {
    if(!file.exists(fp)) {
      mkdirs(dirname(fp))
      dir.create(fp)
    }
  } 
  mkdirs(output.fig.dir)
  mkdirs(output.file.dir)

  ##################################
  # Read in correlated results, selecting just genes as background
  cat("Using the sce file to obtain the background list of genes: '", sce.rds.file, "'")
  #message("\nFiltering thresholds chosen:\nPadj <= ", padj.threshold, "\nand", "\ncorr.est >= ", corr.est.threshold, " OR corr.est <= -", corr.est.threshold)
  #
  message("Reading in sce datafile: '", sce.rds.file, "'")
  #
  data <- readRDS(sce.rds.file)
  # Filter genes for chosen set identified (by correlation analysis)
  # Filtering for only genes
  data_filt <- data[rowData(data)$Type == "gene",]

  # Filter for batch correction method
  message("\nand using the assay/batch corrected data: ",batch.correction.method)
  #
  batch_correction_method <- batch.correction.method
  Gene_Matrix <-   assay(data_filt, batch_correction_method)


  # Using colData as samples sheet, for most up to date QC.
  # Create high vs low methlyation level column (defined methylation level). 
  message("Defining 'high' methylated samples as methylation > ", methylation.threshold, "%")
  #
  colData_filtered <- as.data.frame(colData(data)) %>% 
    mutate(Methylation_level = ifelse(meanMeth > methylation.threshold, "high","low")) %>%
    arrange(desc(Treatment), desc(Methylation_level)) %>%
    dplyr::select(Treatment, Methylation_level) 


  message("Total sample numbers (post-QC) in each treatment and methylation level:")
  print(
    colData_filtered %>% group_by(Treatment, Methylation_level) %>% tally()
    )

  # Should save this table

  #################################################################################################
  # Filter sample IDs for high vs low methylation levels in treatment groups
  message("Filtering into treatment groups, and high vs. low methylation dac.")
  #
  unt <-  colData_filtered %>% filter(Treatment == "unt") %>% row.names(.)
  dac_high  <- colData_filtered %>% filter(Treatment == "dac" & Methylation_level == "high") %>% row.names(.)
  dac_low  <- colData_filtered %>% filter(Treatment == "dac" & Methylation_level == "low") %>% row.names(.)
  aza  <- colData_filtered %>% filter(Treatment == "aza") %>% row.names(.)
  #aza_high  <- colData_filtered %>% filter(Treatment_group == "aza" & Methylation_level == "high") %>% row.names(.)
  #aza_low  <- colData_filtered %>% filter(Treatment_group == "aza" & Methylation_level == "low") %>% row.names(.)


  # Matrix of average values for each treatment group (as identified above)
  message("Creating count matix to heatmap, based oon average count for each treatment group (and high vs low dac).")
  #
  Mean_Matrix <- data.frame(unt =rowMeans(Gene_Matrix[,unt], na.rm = TRUE),
                            dac_high =rowMeans(Gene_Matrix[,dac_high], na.rm = TRUE),
                            dac_low =rowMeans(Gene_Matrix[,dac_low], na.rm = TRUE),
                            aza =rowMeans(Gene_Matrix[,aza], na.rm = TRUE)
  )

  # 2) Row normalised by Mean centering (x -RowMean)
  message("Mean centering data")
  #
  Mean_Matrix_Mean_centred <- Mean_Matrix - rowMeans(Mean_Matrix) 

  message("Preparing heatmapping info...")

  ## Group Info
  Mean_matrix_sample_sheet <- data.frame(Treatment = c("unt", "dac","dac","aza"),
                                         Methylation_level = c("high", "high", "low", "low"))
  row.names(Mean_matrix_sample_sheet) <- colnames(Mean_Matrix)

  # Row annotations
  annotation_row <- data.frame(Type = rowData(data_filt)$Type)
  row.names(annotation_row) <- row.names(data_filt) 


  # Set colours
  Matrix_to_colour_scale <- Mean_Matrix
  paletteLength <- 200
  myColor <- colorRampPalette(c("blue", "light blue", "white", "orange", "red"))(paletteLength)
  myBreaks <- c(seq(min(Matrix_to_colour_scale), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(Matrix_to_colour_scale)/paletteLength, max(Matrix_to_colour_scale), length.out=floor(paletteLength/2)))

  Matrix_to_colour_scale <- Mean_Matrix_Mean_centred
  myBreaks_mean_centred <- c(seq(min(Matrix_to_colour_scale), 0, length.out=ceiling(paletteLength/2) + 1), 
                             seq(max(Matrix_to_colour_scale)/paletteLength, max(Matrix_to_colour_scale), length.out=floor(paletteLength/2)))

  ann_colors = list(
    Treatment_group = c(dac = "red", unt = "green", aza = "blue"),
    Methylation_level = c(low = "black", high = "orchid"),
    Type = c(gene = "green", TE = "blue"))

  ### Run/loop this five times, so the user can see how the clustering algorithm may slightly alter groups upon each run. Seems to repeat every 4-5 runs.
  message("Plotting (and saving results) from 5 repeated heatmaps as clustering lustering algorithm may slightly alter groups upon each run. \nThis seems to stabalised every 4-5 runs.")
  message("\'Stabalisation'/variation seems to depend on the number of genes being plotted AND number of clusters to be defined.")
  cat("\nKmeans (number of groups/clusters) has been set as: ", kmeansk.chosen)

  for(Run in 1:Number.of.plots){
  
    Mean_Matrix_rowMeans_with_kmeans <- 
      pheatmap(Mean_Matrix_Mean_centred,
               kmeans_k = kmeansk.chosen, # decide how many clusters (set to NA if not concerned)
               #treeheight_row = 1,
               color = myColor,
               breaks = myBreaks_mean_centred,
             
               main = paste0(batch.correction.method," - genes Pos. and Neg. Corr. with Methylation in DAC cells \n showing average expression in cell treatment groups (Mean Centred)"),
               show_rownames = T, show_colnames = F,
               cutree_cols = 1,
               clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",  scale = "none", cluster_cols = F,
               annotation_row = annotation_row,
               annotation_col = Mean_matrix_sample_sheet,
               annotation_colors = ann_colors,
               #color = colorRampPalette(c("blue", "light blue", "white", "light orange", "orange", "red"))(30),
             
               annotation_names_row = F,
               fontsize = 8.5)
    #filename = paste0(output.fig.dir,"Heatmap_DAC_Pos.Neg_",Corr.value.tested,"_",batch_correction_method,"_Average_across_treatment_groups_MeanCentred_4Km_Run", Run,".png"))
  
      ggsave(plot=Mean_Matrix_rowMeans_with_kmeans,
            file=paste0(output.fig.dir,"Heatmap_all_genes_",batch_correction_method,"_Average_across_treatment_groups_MeanCentred_",kmeansk.chosen,"Km_Run", Run,".svg"), 
             width=8, height=9)
      ggsave(plot=Mean_Matrix_rowMeans_with_kmeans,
            file=paste0(output.fig.dir,"Heatmap_all_genes_",batch_correction_method,"_Average_across_treatment_groups_MeanCentred_",kmeansk.chosen,"Km_Run", Run,".png"), 
            width=8, height=9)
  


  ### Extract clusters from heatmap ###########################################################################
  # Extract cluster number for each genes and update original correlation table
  ## This bit is more specific to the NHMRC grant

  # Set which pheatmap to extract cluster from
    Pheatmap_to_explore <- Mean_Matrix_rowMeans_with_kmeans
  #Obtain Matrix of clusters
    Gene_Clusters <- as.data.frame(Pheatmap_to_explore$kmeans$cluster) %>% 
      tibble::rownames_to_column( "genes") %>%
      dplyr::rename("cluster_number" = "Pheatmap_to_explore$kmeans$cluster")

  #colnames(Gene_Clusters) <- c("genes","cluster_number")
  #write.csv(Gene_Clusters,
  #          paste0(output.file.dir,"Gene_Clusters_from_Heatmap-Pos.Neg_",Corr.value.tested,"_",batch_correction_method,"_Average_across_treatment_groups_MeanCentred_4.csv"),
  #         row.names = F)

  # Save updated file
    write.csv(Gene_Clusters, 
              file = paste0(output.file.dir,"All_genes_heatmap_cluster_Run",Run,".csv"),
              row.names = F)

  }

}

###########################################################################################################################################################################################################


