##########################################################################
# Purpose: R Script to produce heatmap of all spls selected features
# Output: Heatmap
#
# Date: 15.01.24
# Version: v.0.14.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) Met and Acc heatmaps plot average of sPLS selected windows for each genomic context
## 2) Heatmap column order based on Kmeans clustering of RNA heatmap
## 3) Version showing legend with average of selected sPLS features 
## 4) Produces heatmap with i) hierichical clustering within Kmeans cluster, and ii) columns forced into Tx groups within Kmeans clusters.
# To do: 
##
# Websites:
## 
##########################################################################
# Load packages
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(mixOmics,biomartr,GenomicRanges,GenomicAlignments,
               data.table, ggplot2, RColorBrewer, dplyr, tidyr, tibble, stringr, cowplot, MultiAssayExperiment, pheatmap,
               circlize, ComplexHeatmap, viridis)

#########################################################################################################################################################
# Variables
# Input files
NMT.MAE.file = "./2_results/scNMT/NMT/Autosomal/NMT_MAE.rds"
spls.selected.features.file = "./2_results/scNMT/NMT/Autosomal/pct10_Y100-100_Tx_only/spls_selected_features_metadata.rds"
cell.metadata.file = "./2_results/scNMT/NMT/Autosomal/cell_metadata_updated.tsv"

# Output
output.fig.dir = "./figures/scNMT/NMT/Autosomal/pct10_Y100-100_Tx_only/heatmap/corrected/"

# Additional variables
blocks_met <- c("met_Window3k1k", "met_Promoter", "met_CGI", "met_H3K4me3", "met_H3K27ac")
blocks_acc <- c("acc_Window3k1k", "acc_Promoter", "acc_CGI", "acc_H3K4me3", "acc_H3K27ac")
blocks_rna <- c("rna")
keep_block_order = TRUE
# Comps
comps.to.include <- c(1,2)
#comps.to.include <- c(2)
#########################################################################################################################################################



# Create dir
dir.create(output.fig.dir, showWarnings = FALSE, recursive = TRUE)

# Import data
## Read in input files
cat("\nReading in NMT MAE file: ", NMT.MAE.file)
NMT.MAE <- readRDS(NMT.MAE.file)
cat("\nReading in spls selected features table: ", spls.selected.features.file)
spls.selected.features <- readRDS(spls.selected.features.file)
cat("\nReading in cell metadata file: ", cell.metadata.file)
cell.metadata <- fread(cell.metadata.file) %>% as.data.frame()

# sort cell data into desired order - This needs to then match the matrix used for heatmapping! 
## Aza, DAC, Unt. High to low Meth (Met Window3K1)
cell.metadata$Treatment <- factor(cell.metadata$Treatment, levels = c("AZA", "DAC", "Unt"), ordered = T)
cell.metadata <- cell.metadata %>% 
  arrange(Treatment, desc(met_Window3k1k))



# Met ----------------------------------------------------------------------------------------------------------------------
## Filter for spls selected features. 
## Makes it into a list - each item is an spls block.
Average_spls_selected_features <- TRUE # If true, a single average value is obtained for all selected features. FALSE shows the individual values of selected features.
Average_All_features_entered_into_sPLS_model <- FALSE


cat("\nMet Layer\n")
{ 
  cat("\nCreating heatmap data frame and rowAnnotation.")
  cat("\nFiltering spls for selected features, keeping those in comps:\n", comps.to.include)

  selected_features_list <- NULL
  block <- blocks_met[3]
  for(block in blocks_met){
    # Obtain block
    tmp <- assay(NMT.MAE, block) %>%
      as.data.frame() %>%
      dplyr::select(cell.metadata$sample) # Important to ensure heatmap and heatmap colAnno is correct!!
    # Obtain IDs of selected features for block
    selected_block_feature_IDs <- 
      spls.selected.features %>% 
      filter(block_name == block) %>% 
      filter(comp %in% comps.to.include) %>%
      .[["feature_ID"]]
    # Filter block for selected features -> add to list
    selected_features_list[[block]] <- tmp[row.names(tmp) %in% selected_block_feature_IDs,]
    # Obtain Mean of selected features - optional
    if(Average_spls_selected_features == TRUE){
      selected_features_list[[block]] = colMeans(selected_features_list[[block]], na.rm = T)
      # row.names(selected_features_list[[block]]) <- block
    }
  }
  
  ################################################
  # To use
  ## If wanting legend to contain average of all features entered into sPLS model. Not just the sPLS selected ones....
  if(Average_All_features_entered_into_sPLS_model == TRUE){
    pct_cells_detected = 10 # filter for %detected for average of input. This was the value chosen in the sPLS model.
    keep <- rowSums(!is.na(tmp))/ncol(tmp)*100 >= pct_cells_detected
    tmp_filtered <- tmp[keep,]
  
    tmp_filtered_ColMean <- colMeans(tmp_filtered, na.rm = T)
    tmp_filtered_RowColMean <- sum(tmp_filtered_ColMean) / length(tmp_filtered_ColMean) * 100
  }
  
  ################################################
  
  
  #str(selected_features_list)

  # Filter for selected
  ## Convert list to dataframe for heatmap and simultaneously obtain rowAnnotation df.
  tmp <- NULL
  tmp2 <- data.frame()
  heatmap_df <- NULL
  rowAnnotation_df <- NULL
  block <- blocks_met[3]
  for(block in blocks_met){
    cat("\nUpdating rownames for: ", block)
    tmp <- selected_features_list[[block]]
    # Update rownames.
    if(Average_spls_selected_features == TRUE){
      tmp <- t(as.data.frame(tmp))
      rownames(tmp) <- c(block)
    } else{
      rownames(tmp) <- paste0(block,"_",rownames(tmp)) # This prevents clashing of matching feature_IDs between blocks
    }
    
    heatmap_df <- rbind(heatmap_df, tmp)
  
    cat("\nCreating row annotation df for: ", block)
    tmp2 <- data.frame(block_feature = rownames(tmp)) %>% column_to_rownames(var = "block_feature")
    tmp2 <- tmp2 %>% 
      mutate(block = block) %>% 
      separate(col = block, into = c("molecular_layer", "context"), sep = "_", remove = F) %>%
      mutate(context = ifelse(molecular_layer == "rna", "rna", context))
  
    rowAnnotation_df <- rbind(rowAnnotation_df, tmp2)
  
  
    if(block == blocks_met[length(blocks_met)]){
      cat("\nFinished combining dfs for:\n", blocks_met)
    }
  }
  heatmap_df_met <- heatmap_df
  rowAnnotation_df_met <- rowAnnotation_df
  cat("\nFinished Met layer prep")
}
# Obtain averages
## Can show a feature average per Tx group or if 'FALSE' will show a single average for all samples (combining all Tx groups).

Legened_feature_average_by_Tx <- FALSE
if(Legened_feature_average_by_Tx == "TRUE"){
  
  DAC_met_av <- dplyr::select(as.data.frame(heatmap_df_met), dplyr::contains("DAC")) %>% rowMeans(.)*100
  DAC_met_av <- DAC_met_av %>% as.data.frame(.) %>% rename("Mean_Methylation_DAC" = 1)
  
  AZA_met_av <- dplyr::select(as.data.frame(heatmap_df_met), dplyr::contains("AZA")) %>% rowMeans(.)*100
  AZA_met_av <-  AZA_met_av %>% as.data.frame(.) %>% rename("Mean_Methylation_AZA" = 1)
  
  Unt_met_av <- dplyr::select(as.data.frame(heatmap_df_met), dplyr::contains("Unt")) %>% rowMeans(.)*100
  Unt_met_av <-  Unt_met_av %>% as.data.frame(.) %>% rename("Mean_Methylation_Unt" = 1)
  heatmap_df_met_averages <- cbind(DAC_met_av, AZA_met_av, Unt_met_av)
  
} else if(Legened_feature_average_by_Tx == "FALSE"){
  
  heatmap_df_met_averages <- heatmap_df_met %>% rowMeans(.)*100
  heatmap_df_met_averages <- heatmap_df_met_averages %>% as.data.frame(.) %>% rename("Mean_Methylation" = 1)
  
}

# Acc ----------------------------------------------------------------------------------------------------------------------
## Filter for spls selected features. 
## Makes it into a list - each item is an spls block.

cat("\nAcc Layer\n")
{ 
  cat("\nCreating heatmap data frame and rowAnnotation.")
  cat("\nFiltering spls for selected features, keeping those in comps:\n", comps.to.include)
  
  
  selected_features_list <- NULL

  #block <- blocks_selected[1]
  for(block in blocks_acc){
    # Obtain block
    tmp <- assay(NMT.MAE, block) %>%
      as.data.frame() %>%
      dplyr::select(cell.metadata$sample) # Important to ensure heatmap and heatmap colAnno is correct!!
    # Obtain IDs of selected features for block
    selected_block_feature_IDs <- 
      spls.selected.features %>% 
      filter(block_name == block) %>% 
      filter(comp %in% comps.to.include) %>%
      .[["feature_ID"]]
    # Filter block for selected features -> add to list
    selected_features_list[[block]] <- tmp[row.names(tmp) %in% selected_block_feature_IDs,]
    # Obtain Mean of selected features - optional
    if(Average_spls_selected_features == TRUE){
      selected_features_list[[block]] = colMeans(selected_features_list[[block]], na.rm = T)
      # row.names(selected_features_list[[block]]) <- block
    }
  }

  #selected_features_list

  ## Convert list to dataframe for heatmap and simultaneously obtain rowAnnotation df.
  tmp <- NULL
  tmp2 <- data.frame()
  heatmap_df <- NULL
  rowAnnotation_df <- NULL
  block <- blocks_acc[1]
  for(block in blocks_acc){
    cat("\nUpdating rownames for: ", block)
    tmp <- selected_features_list[[block]]
    # Update rownames.
    if(Average_spls_selected_features == TRUE){
      tmp <- t(as.data.frame(tmp))
      rownames(tmp) <- c(block)
    } else{
      rownames(tmp) <- paste0(block,"_",rownames(tmp)) # This prevents clashing of matching feature_IDs between blocks
    }
    
    heatmap_df <- rbind(heatmap_df, tmp)
  
    cat("\nCreating row annotation df for: ", block)
    tmp2 <- data.frame(block_feature = rownames(tmp)) %>% column_to_rownames(var = "block_feature")
    tmp2 <- tmp2 %>% 
      mutate(block = block) %>% 
      separate(col = block, into = c("molecular_layer", "context"), sep = "_", remove = F) %>%
      mutate(context = ifelse(molecular_layer == "rna", "rna", context))
  
    rowAnnotation_df <- rbind(rowAnnotation_df, tmp2)
  
  
    if(block == blocks_acc[length(blocks_acc)]){
      cat("\nFinished combining dfs for:\n", blocks_acc)
    }
  }
  heatmap_df_acc <- heatmap_df
  rowAnnotation_df_acc <- rowAnnotation_df
  cat("\nFinished Acc layer prep")
}
# Obtain averages
## Can show a feature average per Tx group or if 'FALSE' will show a single average for all samples (combining all Tx groups).

Legened_feature_average_by_Tx <- FALSE
if(Legened_feature_average_by_Tx == TRUE){
  
  DAC_acc_av <- dplyr::select(as.data.frame(heatmap_df_acc), dplyr::contains("DAC")) %>% rowMeans(.)*100
  DAC_acc_av <- DAC_acc_av %>% as.data.frame(.) %>% rename("Mean_Accessibility_DAC" = 1)
  
  AZA_acc_av <- dplyr::select(as.data.frame(heatmap_df_acc), dplyr::contains("AZA")) %>% rowMeans(.)*100
  AZA_acc_av <-  AZA_acc_av %>% as.data.frame(.) %>% rename("Mean_Accessibility_AZA" = 1)
  
  Unt_acc_av <- dplyr::select(as.data.frame(heatmap_df_acc), dplyr::contains("Unt")) %>% rowMeans(.)*100
  Unt_acc_av <-  Unt_acc_av %>% as.data.frame(.) %>% rename("Mean_Accessibility_Unt" = 1)
  heatmap_df_acc_averages <- cbind(DAC_acc_av, AZA_acc_av, Unt_acc_av)
  
} else if(Legened_feature_average_by_Tx == FALSE){
  
  heatmap_df_acc_averages <- heatmap_df_acc %>% rowMeans(.)*100
  heatmap_df_acc_averages <- heatmap_df_acc_averages %>% as.data.frame(.) %>% rename("Mean_Accessibility" = 1)

}
# RNA ----------------------------------------------------------------------------------------------------------------------
## Filter for spls selected features. 
## Makes it into a list - each item is an spls block.

cat("\nRNA Layer\n")
{ 
  cat("\nCreating heatmap data frame and rowAnnotation.")
  cat("\nFiltering spls for selected features, keeping those in comps:\n", comps.to.include)
  
  selected_features_list <- NULL
  block <- blocks_rna[1]
  for(block in blocks_rna){
    # Obtain block
    tmp <- assay(NMT.MAE, block) %>%
      as.matrix() %>%
      as.data.frame() %>%
      dplyr::select(cell.metadata$sample) # Important to ensure heatmap and heatmap colAnno is correct!!
    # Obtain IDs of selected features for block
    selected_block_feature_IDs <- 
      spls.selected.features %>% 
      filter(block_name == block) %>% 
      filter(comp %in% comps.to.include) %>%
      .[["feature_ID"]]
    # Filter block for selected features -> add to list
    selected_features_list[[block]] <- tmp[row.names(tmp) %in% selected_block_feature_IDs,]
 }
  ## Convert list to dataframe for heatmap and simultaneously obtain rowAnnotation df.
  tmp <- NULL
  tmp2 <- data.frame()
  heatmap_df <- NULL
  rowAnnotation_df <- NULL
  block <- blocks_rna[1]
  for(block in blocks_rna){
    cat("\nUpdating rownames for: ", block)
    tmp <- selected_features_list[[block]]
    rownames(tmp) <- paste0(block,"_",rownames(tmp)) # This prevents clashing of matching feature_IDs between blocks
  
    heatmap_df <- rbind(heatmap_df, tmp)
  
    cat("\nCreating row annotation df for: ", block)
    tmp2 <- data.frame(block_feature = rownames(tmp)) %>% column_to_rownames(var = "block_feature")
    tmp2 <- tmp2 %>% 
      mutate(block = block) %>% 
      separate(col = block, into = c("molecular_layer", "context"), sep = "_", remove = F) %>%
      mutate(context = ifelse(molecular_layer == "rna", "rna", context))
  
    rowAnnotation_df <- rbind(rowAnnotation_df, tmp2)
  
  if(block == blocks_rna[length(blocks_rna)]){
      cat("\nFinished combining dfs for:\n", blocks_rna)
      }
  }
  heatmap_df_rna <- heatmap_df
  rowAnnotation_df_rna <- rowAnnotation_df
  cat("\nFinished RNA layer prep")
}
# Heatmap 1 - Acc #########################################################################################################################################
##################
convert_to_mean_centred = FALSE
convert_to_z_score = TRUE # If FALSE, the 'raw' values will be used (no row normalisation)
#Legened_feature_average_by_Tx = FALSE

cat("\nStarting ComplexHeatmap for 'met' blocks:\n",blocks_acc)
{
  df_to_heatmap <- heatmap_df_acc
  blocks_selected <- blocks_acc
  
  if(convert_to_z_score|convert_to_mean_centred == TRUE){
    df_to_heatmap_raw <- df_to_heatmap
    df_to_heatmap_raw_Mean_centred <- df_to_heatmap_raw - rowMeans(df_to_heatmap_raw, na.rm = T)
    df_to_heatmap_Zscore <- t(scale(t(df_to_heatmap_raw_Mean_centred)))
    if(convert_to_z_score == TRUE) {
      cat("\n'df_to_heatmap' has been converted to z-score.\n")
      df_to_heatmap <- df_to_heatmap_Zscore
    } else if(convert_to_mean_centred == TRUE){
      cat("\n'df_to_heatmap' has been converted to mean centred.\n")
      df_to_heatmap <- df_to_heatmap_raw_Mean_centred
    }
    cat("\n'The raw df can be accessed via 'df_to_heatmap_raw' and \nz-score can also be accessed via 'df_to_heatmap_Zscore'.\n")
  }


  colAnnotation_df <- cell.metadata %>% 
    dplyr::select(sample,Treatment) %>% # , all_of( blocks_selected[!grepl("rna",blocks_selected)]) 
    column_to_rownames("sample")

  # Make sure using correct col dataframes and Matrix for heatmap!
  ## Column info
  #  colAnnotation_df
  ## Left annotation (spls blocks)
  rowAnnotation_df_left <- rowAnnotation_df_acc %>% dplyr::select(-block)
  rowAnnotation_df_left$context <- factor(rowAnnotation_df_left$context, levels = c("Window3k1k", "Promoter", "CGI", "H3K4me3","H3K27ac"))
  ## Right annotation (average values of rows)
  rowAnnotation_df_right <- heatmap_df_acc_averages # %>% as.data.frame(.) %>% rename("Mean_Accesibility" = 1)
  #df_to_heatmap

  # Create colour panels/df
  ## Matrix colour
  min_val = min(df_to_heatmap, na.rm =T)
  max_val = max(df_to_heatmap, na.rm = T)


  ## Column Anno colour
  anno_colors_col = list(
    Treatment= c(DAC = "#00FFFF", Unt = "orange", AZA = "orchid")
  )
#  for(blocks in blocks_selected){
#    anno_colors_col[[blocks]] = colour_fun_ColAnno 
#  }

  ## Row Anno colours
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  ### Block colours
  block_colours <- col_vector[1:length(blocks_selected)]
  names(block_colours) <- blocks_selected
  #### context colours
  context_colours <- col_vector[30:(29+length(unique(rowAnnotation_df_left$context)))]
  names(context_colours) <- unique(rowAnnotation_df_left$context)
  #### Left row (layers)
  anno_colors_row_left = list(
    block = block_colours,
    molecular_layer = c(met = "brown", acc = "blue", rna = "green"),
    context = context_colours)
  #### right row (average levels)
  
  if(Legened_feature_average_by_Tx == TRUE){
    anno_colors_row_right = list(
      Mean_Accessibility_DAC = colorRamp2(c(0,65, 100), c(inferno(300)[1],inferno(300)[200], inferno(300)[296])),
      Mean_Accessibility_AZA = colorRamp2(c(0,65, 100), c(inferno(300)[1],inferno(300)[200], inferno(300)[296])),
      Mean_Accessibility_Unt = colorRamp2(c(0,65, 100), c(inferno(300)[1],inferno(300)[200], inferno(300)[296])))
    
    anno_colors_row_right = list(
      Mean_Accessibility_DAC = colorRamp2(c(0,30, 60, 100), c(inferno(300)[1],inferno(300)[170],inferno(300)[250], inferno(300)[296])),
      Mean_Accessibility_AZA = colorRamp2(c(0,30, 60, 100), c(inferno(300)[1],inferno(300)[170],inferno(300)[250], inferno(300)[296])),
      Mean_Accessibility_Unt = colorRamp2(c(0,30, 60, 100), c(inferno(300)[1],inferno(300)[170],inferno(300)[250], inferno(300)[296])))
    
    
  } else if(Legened_feature_average_by_Tx == FALSE){
    anno_colors_row_right = list(
    Mean_Accessibility = colorRamp2(c(0,65, 100), c(inferno(300)[1],inferno(300)[200], inferno(300)[296])))
  }
  
  # Convert to ComplexHeatmap Annotation info
  ha_col <- ComplexHeatmap::HeatmapAnnotation(df = colAnnotation_df, which = "col", col = anno_colors_col,
                                              show_legend = c(TRUE, TRUE, rep(FALSE,5)))
  ha_row_left <- ComplexHeatmap::HeatmapAnnotation(df = rowAnnotation_df_left, which = "row", col = anno_colors_row_left, 
                                              annotation_name_rot = 90,
                                              show_annotation_name = T, annotation_name_side = "bottom")
  ha_row_right <- ComplexHeatmap::HeatmapAnnotation(df = rowAnnotation_df_right, which = "row", col = anno_colors_row_right, 
                                                   annotation_name_rot = 90,
                                                   show_annotation_name = T, annotation_name_side = "bottom", annotation_label = "Mean_Level")


  colour_fun_Matrix = colorRamp2(c(0, 0.3, 0.4, 0.6, 0.8, 1), c("#BFFFC3", "#4FC857","yellow", "orange", "red", "darkred"))
  # Put into heatmap
  if(convert_to_z_score|convert_to_mean_centred == TRUE){
    colour_fun_Matrix = colorRamp2(c(min_val,0-1e-8, 0,0+1e-08,max_val/4, max_val), c("darkblue", "lightblue", "white","#FFDEC8", "#FF9D9D", "red"))
  }

  
  # Plot Heatmap
  CH_acc <- 
    ComplexHeatmap::Heatmap(matrix = as.matrix(df_to_heatmap),
                            name = "Accessibility level",
                            # col = viridis(500), 
                            col = colour_fun_Matrix,
                            show_column_names = F, show_column_dend = T,
                            show_row_names = F, row_names_side = "left",
                            show_row_dend = F,
                        
                            cluster_rows = T, cluster_columns = T,
                            row_split = rowAnnotation_df_left$context, cluster_row_slices = F,
                            #top_annotation = ha_col,
                            left_annotation = ha_row_left,
                            right_annotation = ha_row_right,
                            heatmap_height = unit(8, "cm"),
                            heatmap_width = unit(6, "cm"),
                            row_title_rot = 0)
  CH_acc
  # Store values uniquely so they can be reordered later
  df_to_heatmap_acc <- as.matrix(df_to_heatmap)
  colour_fun_Matrix_acc <- colour_fun_Matrix
  rowAnnotation_df_left_acc <- rowAnnotation_df_left
  ha_row_left_acc <- ha_row_left
  right_annotation_acc <- ha_row_right
  rowAnnotation_df_right_acc <- rowAnnotation_df_right
}

# Heatmap 2 - Met #########################################################################################################################################
##################
convert_to_mean_centred = FALSE
convert_to_z_score = TRUE # If FALSE, the 'raw' values will be used (no row normalisation)

cat("\nStarting ComplexHeatmap for 'met' blocks:\n",blocks_met)
{
  df_to_heatmap <- heatmap_df_met
  blocks_selected <- blocks_met
  
  if(convert_to_z_score|convert_to_mean_centred == TRUE){
    df_to_heatmap_raw <- df_to_heatmap
    df_to_heatmap_raw_Mean_centred <- df_to_heatmap_raw - rowMeans(df_to_heatmap_raw, na.rm = T)
    df_to_heatmap_Zscore <- t(scale(t(df_to_heatmap_raw_Mean_centred)))
    if(convert_to_z_score == TRUE) {
      cat("\n'df_to_heatmap' has been converted to z-score.\n")
      df_to_heatmap <- df_to_heatmap_Zscore
    } else if(convert_to_mean_centred == TRUE){
      cat("\n'df_to_heatmap' has been converted to mean centred.\n")
      df_to_heatmap <- df_to_heatmap_raw_Mean_centred
    }
    cat("\n'The raw df can be accessed via 'df_to_heatmap_raw' and \nz-score can also be accessed via 'df_to_heatmap_Zscore'.\n")
  }

  colAnnotation_df <- cell.metadata %>% 
    dplyr::select(sample,Treatment) %>% # all_of( blocks_selected[!grepl("rna",blocks_selected)]) 
    column_to_rownames("sample")

  # Make sure using correct col dataframes and Matrix for heatmap!
  colAnnotation_df
  ## Left annotation (spls blocks)
  rowAnnotation_df_left <- rowAnnotation_df_met %>% dplyr::select(-block)
  rowAnnotation_df_left$context <- factor(rowAnnotation_df_left$context, levels = c("Window3k1k", "Promoter", "CGI", "H3K4me3","H3K27ac"))
  ## Right annotation (average values of rows)
  rowAnnotation_df_right <- heatmap_df_met_averages # %>% as.data.frame(.) %>% rename("Mean_Methylation" = 1)

  #df_to_heatmap

  # Create colour panels/df
  ## Matrix colour
  min_val = min(df_to_heatmap, na.rm =T)
  max_val = max(df_to_heatmap, na.rm = T)

  ## Column Anno colour
#  colour_fun_ColAnno = colorRamp2(c(0, 50,50+1e-6,100), c("white", "blue","lightgreen", "darkgreen"))

  anno_colors_col = list(
    Treatment= c(DAC = "#00FFFF", Unt = "orange", AZA = "orchid")
    )
#  for(blocks in blocks_selected){
#    anno_colors_col[[blocks]] = colour_fun_ColAnno 
#  }

  ## Row Anno colours
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  ### Block colours
  block_colours <- col_vector[1:length(blocks_selected)]
  names(block_colours) <- blocks_selected
  #### context colours
  context_colours <- col_vector[30:(29+length(unique(rowAnnotation_df_left$context)))]
  names(context_colours) <- unique(rowAnnotation_df_left$context)

  anno_colors_row_left = list(
 #   block = block_colours,
    molecular_layer = c(met = "brown", acc = "blue", rna = "green"),
    context = context_colours)
  
  
  if(Legened_feature_average_by_Tx == TRUE){
    
    anno_colors_row_right = list(
      Mean_Methylation_DAC = colorRamp2(c(0,30, 60, 100), c(inferno(300)[1],inferno(300)[170],inferno(300)[250], inferno(300)[296])),
      Mean_Methylation_AZA = colorRamp2(c(0,30, 60, 100), c(inferno(300)[1],inferno(300)[170],inferno(300)[250], inferno(300)[296])),
      Mean_Methylation_Unt = colorRamp2(c(0,30, 60, 100), c(inferno(300)[1],inferno(300)[170],inferno(300)[250], inferno(300)[296])))
    
    
  } else if(Legened_feature_average_by_Tx == FALSE){
  
  anno_colors_row_right = list(
    Mean_Methylation = colorRamp2(c(0,65, 100), c(inferno(300)[1],inferno(300)[200], inferno(300)[296])))
  }
    #Mean_Methylation = colorRamp2(c(0,50,70, 100), c(inferno(300)[1], inferno(300)[200], inferno(300)[250], inferno(300)[296]))
    #Mean_Methylation = colorRamp2(c(0, 70,70+1e-6,100), c("white", "darkblue","#FFD2D2", "#720000"))



  # Convert to ComplexHeatmap Annotation info
  ha_col <- ComplexHeatmap::HeatmapAnnotation(df = colAnnotation_df, which = "col", col = anno_colors_col,
                                              show_legend = c(TRUE, TRUE, rep(FALSE,5)))
  ha_row_left <- ComplexHeatmap::HeatmapAnnotation(df = rowAnnotation_df_left, which = "row", col = anno_colors_row_left, 
                                              show_annotation_name = F, annotation_name_side = "top",
                                              annotation_name_rot = 90,
                                              show_legend = c(TRUE, TRUE, rep(FALSE,4),TRUE,rep(FALSE,4)))
  ha_row_right <- ComplexHeatmap::HeatmapAnnotation(df = rowAnnotation_df_right, which = "row", col = anno_colors_row_right, 
                                                    annotation_name_rot = 90,
                                                    show_annotation_name = F, annotation_name_side = "top") #, annotation_label = "Mean_Level")


  colour_fun_Matrix = colorRamp2(c(0, 0.3, 0.7, 1), c("lightgreen","darkgreen", "black", "red"))
  if(convert_to_z_score|convert_to_mean_centred == TRUE){
    colour_fun_Matrix = colorRamp2(c(min_val,0-1e-8, 0,0+1e-08,max_val/4, max_val), c("darkblue", "lightblue", "white","#FFDEC8", "#FF9D9D", "red"))
  } 

  # Put into heatmap
  CH_met <- 
  ComplexHeatmap::Heatmap(matrix = as.matrix(df_to_heatmap),
                          name = "DNA methylation level",
                          col = colour_fun_Matrix,
                          show_column_names = F, show_column_dend = T,
                          show_row_names = F, row_names_side = "left",
                          show_row_dend = F,
                        
                          cluster_rows = T, cluster_columns = T,
                          row_split = rowAnnotation_df_left$context, cluster_row_slices = F,
                          # top_annotation = ha_col,
                          left_annotation = ha_row_left,
                          right_annotation = ha_row_right,
                          heatmap_height = unit(8, "cm"),
                          heatmap_width = unit(6, "cm"),
                          row_title_rot = 0
                          )
  CH_met
  # Store values uniquely so they can be reordered later
  df_to_heatmap_met <- as.matrix(df_to_heatmap)
  colour_fun_Matrix_met <- colour_fun_Matrix
  rowAnnotation_df_left_met <- rowAnnotation_df_left
  ha_row_left_met <- ha_row_left
  right_annotation_met <- ha_row_right
  rowAnnotation_df_right_met <- rowAnnotation_df_right
}

# Heatmap 3 - RNA #########################################################################################################################################
##################
##################
convert_to_z_score = TRUE # If FALSE, the 'raw' values will be used (no row normalisation)

cat("\nStarting ComplexHeatmap for 'met' blocks:\n",blocks_acc)
{
  df_to_heatmap <- heatmap_df_rna
  blocks_selected <- blocks_rna

  if(convert_to_z_score == TRUE){
    df_to_heatmap_raw <- df_to_heatmap
    df_to_heatmap_raw_Mean_centred <- df_to_heatmap_raw - rowMeans(df_to_heatmap_raw, na.rm = T)
    df_to_heatmap_Zscore <- t(scale(t(df_to_heatmap_raw_Mean_centred)))
    df_to_heatmap <- df_to_heatmap_Zscore
    cat("\n'df_to_heatmap' has been converted to z-score.\n")
    cat("\n'The raw df can be accessed via 'df_to_heatmap_raw' and \nz-score can also be accessed via 'df_to_heatmap_Zscore'.\n")
  }


  colAnnotation_df <- cell.metadata %>% 
    dplyr::select(sample,Treatment) %>% #all_of( blocks_selected[!grepl("rna",blocks_selected)]) 
    column_to_rownames("sample")

  # Make sure using correct col dataframes and Matrix for heatmap!
  colAnnotation_df
  # Left row annotation
  rowAnnotation_df_left <- rowAnnotation_df_rna %>% dplyr::select(-block)
  # right row annotation
  rowAnnotation_df_right <- rowAnnotation_df_left %>% 
    mutate(Type = ifelse(grepl(pattern = "ENSG", x = rownames(rowAnnotation_df_left)), "Gene", "TE")) %>%
    dplyr::select(Type)

  # Create colour panels/df
  ## Matrix colour
  min_val = min(df_to_heatmap, na.rm =T)
  max_val = max(df_to_heatmap, na.rm = T)
  colour_fun_Matrix = colorRamp2(c(min_val, 0, max_val), c("darkblue", "white", "darkorange"))
  #colour_fun_Matrix = colorRamp2(c(0, 0.5, 1), c("lightgreen", "black", "red"))
  
  ## Column Anno colour
  colour_fun_ColAnno = colorRamp2(c(0, 50,50+1e-6,100), c("white", "blue","lightgreen", "darkgreen"))

  anno_colors_col = list(
    Treatment= c(DAC = "#00FFFF", Unt = "orange", AZA = "orchid")
  )
  for(blocks in blocks_selected){
    anno_colors_col[[blocks]] = colour_fun_ColAnno 
  }

  ## Row Anno colours
  n <- 60
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  ### Block colours
  block_colours <- col_vector[1:length(blocks_selected)]
  names(block_colours) <- blocks_selected
  #### context colours
  context_colours <- col_vector[30:(29+length(unique(rowAnnotation_df$context)))]
  names(context_colours) <- unique(rowAnnotation_df$context)

  anno_colors_row_left = list(
    block = block_colours,
    molecular_layer = c(met = "brown", acc = "blue", rna = "#00B80E"),
    context = context_colours)
  
  
  anno_colors_row_right = list(
    Type = c(Gene = "darkblue", TE = "yellow"))

  # Convert to ComplexHeatmap Annotation info
  ha_col <- ComplexHeatmap::HeatmapAnnotation(df = colAnnotation_df, which = "col", col = anno_colors_col, annotation_name_side = "left")
  ha_row_left <- ComplexHeatmap::HeatmapAnnotation(df = rowAnnotation_df_left, which = "row", col = anno_colors_row_left, 
                                              show_annotation_name = F, annotation_name_side = "top")
  ha_row_right <- ComplexHeatmap::HeatmapAnnotation(df = rowAnnotation_df_right, which = "row", col = anno_colors_row_right, 
                                              show_annotation_name = T, annotation_name_side = "top")


  colour_fun_Matrix = colorRamp2(c(min_val,0-1e-6, 0,0+1e-06, max_val), c("darkblue", "lightblue", "white","#FFDEC8", "darkred"))
  # Heatmap
  CH_rna_splitting_kmeans <-
    ComplexHeatmap::Heatmap(matrix = as.matrix(df_to_heatmap),
                            name = "RNA z-score",
                            col = colour_fun_Matrix,
                            show_column_names = F, show_column_dend = T,
                            show_row_names = F, show_row_dend = F,
                        
                            cluster_rows = T, cluster_columns = T,
                            # row_split = rowAnnotation_df$context, 
                            cluster_row_slices = T,
                            top_annotation = ha_col,
                            left_annotation = ha_row_left,
                            right_annotation = ha_row_right,
                            row_km = 3)
  
  # Store values uniquely so they can be reordered later
  df_to_heatmap_rna<- as.matrix(df_to_heatmap)
  colour_fun_Matrix_rna <- colour_fun_Matrix
  rowAnnotation_df_left_rna <- rowAnnotation_df_left
  ha_row_left_rna <- ha_row_left
  right_annotation_rna <- ha_row_right
}


# Combine Heatmap 1-3 #####################################################################################################################################
######################
cat("Combining heatmaps")
ht_list_rna_kmeans_split <- CH_met %v% CH_acc %v% CH_rna_splitting_kmeans

if(length(comps.to.include) > 1){
  comps.to.include_name <- paste(comps.to.include,collapse='-')
} else{
  comps.to.include_name <- comps.to.include
}

set.seed(123)

output.file.name <- paste0(output.fig.dir, "Heatmap_NMT_C", comps.to.include_name,"_RNA_rowKmeans3_colKmeans4.pdf")
cat("\nSaving final heatmap into: \n", output.file.name)
pdf(output.file.name, width = 12, height = 15)
rowKmeans3_colKmeans4 <- draw(ht_list_rna_kmeans_split[c(3,1,2)], main_heatmap = "RNA z-score", column_km = 4,column_km_repeats = 101)
dev.off()


######################
# Obtain IDs from clusters for samples (columns) and RNA (rows).
## Samples (columns) ------------------------------------------------------------------------------------------------------------------------------
cat("\nExtracting Kmeans clusters from...\n ") 
cat("\nColumns (Samples)\n ")
Sample_Cluster <- NULL
df <- NULL
df_tmp <-NULL
for(cluster in names(column_order(rowKmeans3_colKmeans4))){
  Sample_Cluster[[cluster]] <- colnames(df_to_heatmap)[column_order(rowKmeans3_colKmeans4)[[cluster]]]
  sample_IDs <- colnames(df_to_heatmap)[column_order(rowKmeans3_colKmeans4)[[cluster]]]
  df_tmp <- data.frame("sample_ID" = sample_IDs,
                       "cluster" = rep(cluster, times = length(sample_IDs)))
  df <-  rbind(df, df_tmp)
}
# Update df

df_upd <-
  df %>% mutate(Treatment = ifelse(grepl(pattern = "AZA", x = sample_ID), "AZA", 
                                 ifelse(grepl(pattern = "DAC", x = sample_ID), "DAC",  
                                        ifelse(grepl(pattern = "Unt", x = sample_ID), "Unt", "NA"))),
              Batch =     ifelse(grepl(pattern = "b620", x = sample_ID), "b620", 
                                 ifelse(grepl(pattern = "b320", x = sample_ID), "b320",  
                                        ifelse(grepl(pattern = "b322", x = sample_ID), "b322", "NA")))
              
              )
df_upd_cols <- df_upd
# Save
output.file.name <- paste0(output.fig.dir, "Heatmap_NMT_C", comps.to.include_name,"_RNA_rowKmeans3_colKmeans4-ClusterIDs_SAMPLES.csv")
cat("\nSaving into: ", output.file.name)
write.csv(df_upd, file = output.file.name, row.names = F)

## Genes/TEs (rows) ------------------------------------------------------------------------------------------------------------------------------
cat("\nRows (RNA)\n ")
Sample_Cluster <- NULL
df <- NULL
df_tmp <-NULL
for(cluster in names(row_order(rowKmeans3_colKmeans4)[['RNA z-score']])){
  Sample_Cluster[[cluster]] <- rownames(df_to_heatmap)[row_order(rowKmeans3_colKmeans4)[['RNA z-score']][[cluster]]]
  sample_IDs <- rownames(df_to_heatmap)[row_order(rowKmeans3_colKmeans4)[['RNA z-score']][[cluster]]]
  df_tmp <- data.frame("row_ID" = sample_IDs,
                       "cluster" = rep(cluster, times = length(sample_IDs)))
  df <-  rbind(df, df_tmp)
}
df_upd <- 
  df %>% separate(col = row_ID, into = c("temp", "feature_ID"), sep = "_", remove = F) %>%
  dplyr::select(row_ID, cluster, feature_ID) # %>% mutate(Type = ifelse(grepl(pattern = "ENSG", x = ENSG.TE), "Gene", "TE"))

# Add hgnc symbol info
rna_metadata <- spls.selected.features %>% filter(block_name == "rna") 
df_upd <-  left_join(df_upd, rna_metadata, by = c("feature_ID" = "feature_ID"))
df_upd_rows <- df_upd
# Save
output.file.name <- paste0(output.fig.dir, "Heatmap_NMT_C", comps.to.include_name,"_RNA_rowKmeans3_colKmeans4-ClusterIDs_GENES.csv")
cat("\nSaving into: ", output.file.name)
write.csv(df_upd, file = output.file.name, row.names = F)


cat("\nHuzzah! It's finally finished. If this has worked, you're a genius! (If not, make it work to become one...)")


############## Re-draw to column order matches desired order.
# (Re)Order cols (samples) by cluster order
## Order cluster df.
df_upd_cols$cluster <- factor(df_upd_cols$cluster, levels = c("4","3","2","1"))
df_upd_cols <- df_upd_cols %>% arrange(cluster, Treatment) 
colAnnotation_df2 <- df_upd_cols %>% column_to_rownames("sample_ID") %>% select(Treatment)
## Update heatmap anno
ha_col <- ComplexHeatmap::HeatmapAnnotation(df = colAnnotation_df2, which = "col", col = anno_colors_col, annotation_name_side = "left")
## Re-order rna matrix
df_to_heatmap_rna2 <- as.data.frame(df_to_heatmap_rna) %>% dplyr::select(df_upd_cols$sample_ID) 

# (Re)Order rows (genes) by cluster order
## Order cluster df.
df_upd_rows$cluster <- factor(df_upd_rows$cluster, levels = c("3","2","1"))
df_upd_rows <- df_upd_rows %>% arrange(cluster) %>% 
  dplyr::distinct(row_ID, .keep_all = TRUE) ## Need to remove duplicate rows created from multiple hgnc symbol per ensembl ID
rowAnnotation_df2 <- df_upd_rows %>% column_to_rownames("row_ID") %>% select(Type) %>%
  mutate(Type = gsub(x = Type, pattern = "gene", replacement = "Gene")) # correct "gene" to "Gene"
## Re-order rna matrix to match cluster order # https://stackoverflow.com/questions/55449896/sorting-rownames-of-a-data-frame-based-on-a-separate-vector
df_to_heatmap_rna2 <- df_to_heatmap_rna2[order(match(rownames(df_to_heatmap_rna2), df_upd_rows$row_ID)), , drop = FALSE]


# Update annos
## Update right anno
ha_row_right2 <- ComplexHeatmap::HeatmapAnnotation(df = rowAnnotation_df2, which = "row", col = anno_colors_row_right, annotation_name_side = "top")
## Update left anno
rowAnnotation_df_left2 <- data.frame(row.name = df_upd_rows$row_ID, molecular_layer = "rna", context = "rna") %>% 
  column_to_rownames("row.name")
ha_row_left2 <- ComplexHeatmap::HeatmapAnnotation(df = rowAnnotation_df_left2, which = "row", col = anno_colors_row_left, 
                                                 show_annotation_name = F, annotation_name_side = "top")

CH_rna_reordered <- 
  ComplexHeatmap::Heatmap(matrix = as.matrix(df_to_heatmap_rna2),
                        name = "RNA z-score",
                        col = colour_fun_Matrix,
                        show_column_names = F, show_column_dend = T,
                        show_row_names = F, show_row_dend = F,
                        
                        cluster_rows = F, cluster_columns = F,
                        # row_split = rowAnnotation_df$context, 
                        cluster_row_slices = T,
                        top_annotation = ha_col,
                        left_annotation = ha_row_left2,
                        right_annotation = ha_row_right2,
                        row_split = df_upd_rows$cluster,
                        column_split = df_upd_cols$cluster)


# Re-order Met heatmap #########################################
## by rna sample kmeans clusters
df_to_heatmap_met.reordered <- as.data.frame(df_to_heatmap_met) %>% dplyr::select(df_upd_cols$sample_ID)  ## by rna sample kmeans clusters

anno_colors_row_right2 = list(
  Mean_Methylation = colorRamp2(c(0,65, 100), c(inferno(300)[1],inferno(300)[200], inferno(300)[296])))

anno_colors_row_right2 = list(
  Mean_Methylation = colorRamp2(c(0,15, 30, 50, 100), c(inferno(100)[1],inferno(100)[15],inferno(100)[50], inferno(100)[75], inferno(100)[100])))

ha_row_right2 <- ComplexHeatmap::HeatmapAnnotation(df = rowAnnotation_df_right_met, which = "row", col = anno_colors_row_right2, 
                                                  annotation_name_rot = 90,
                                                  show_annotation_name = F, annotation_name_side = "top") #, annotation_label = "Mean_Level")


CH_met_reordered <-
  ComplexHeatmap::Heatmap(matrix = as.matrix(df_to_heatmap_met.reordered),
                        #column_order = df_upd_cols$sample_ID, ## by rna sample kmeans clusters
                        
                        name = "DNA methylation level",
                        col = colour_fun_Matrix_met,
                        show_column_names = F, show_column_dend = T,
                        show_row_names = F, row_names_side = "left",
                        show_row_dend = F,
                        
                        cluster_rows = T, cluster_columns = F,
                        row_split = rowAnnotation_df_left_met$context, cluster_row_slices = F,
                        # top_annotation = ha_col,
                        left_annotation = ha_row_left_met,
                        right_annotation = ha_row_right2,
                        heatmap_height = unit(8, "cm"),
                        heatmap_width = unit(6, "cm"),
                        row_title_rot = 0)


# Re-order Acc heatmap ##############################################
## by rna sample kmeans clusters
df_to_heatmap_acc.reordered <- as.data.frame(df_to_heatmap_acc) %>% dplyr::select(df_upd_cols$sample_ID)  ## by rna sample kmeans clusters


anno_colors_row_right2 = list(
  Mean_Accessibility = colorRamp2(c(0,15, 30, 50, 100), c(inferno(100)[1],inferno(100)[15],inferno(100)[50], inferno(100)[75], inferno(100)[100])))

ha_row_right2_acc <- ComplexHeatmap::HeatmapAnnotation(df = rowAnnotation_df_right_acc, which = "row", col = anno_colors_row_right2, 
                                                   annotation_name_rot = 90,
                                                   show_annotation_name = F, annotation_name_side = "top") #, annotation_label = "Mean_Level")

CH_acc_reordered <-
  ComplexHeatmap::Heatmap(matrix = as.matrix(df_to_heatmap_acc.reordered),
                          #column_order = df_upd_cols$sample_ID, ## by rna sample kmeans clusters
                          
                          
                          name = "Accessibility level",
                          # col = viridis(500), 
                          col = colour_fun_Matrix_acc,
                          show_column_names = F, show_column_dend = T,
                          show_row_names = F, row_names_side = "left",
                          show_row_dend = F,
                          
                          cluster_rows = T, cluster_columns = T,
                          row_split = rowAnnotation_df_left_acc$context, cluster_row_slices = F,
                          #top_annotation = ha_col,
                          left_annotation = ha_row_left_acc,
                          right_annotation = ha_row_right2_acc,
                          heatmap_height = unit(8, "cm"),
                          heatmap_width = unit(6, "cm"),
                          row_title_rot = 0)

# Combine updated ##################################
ht_list_updated <- CH_rna_reordered %v% CH_met_reordered %v% CH_acc_reordered

output.file.name <- paste0(output.fig.dir, "Heatmap_NMT_C", comps.to.include_name,"_RNA_rowKmeans3_colKmeans4_Tx_ordered.pdf")
cat("\nSaving final reordered heatmap into: \n", output.file.name)
pdf(output.file.name, width = 12, height = 15)
  draw(ht_list_updated)
dev.off()
############################################

