##########################################################################
# Purpose: R Script to filter sce samples based on prior seqmonk QC values (requires 'Passed_seqmonk_QC' column in metadata) and produce relevant QC plots.
# Output: i) filtered sce .rds files and ii) several QC plots pre (QC1) and post (QC2) filtered
## sce_qc2.rds is the QC (sample) + feature filtered file, while sce_qc1.rds is simply the imported sce.rds + initial QC metrics (no feature filtering).
#
# Date: 02.Feb.23
# Version: v.0.2.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
##  1) Contains functions i) perform_QC_seqmonk_based and ii) perform_QC_scRNAseq. Although both rely on seqmonk QC filtering table in the meta data... *******************
##  2) Filtering of lowly expressed genes by min.counts in x cells (default) or x% of cells.
##  3) Corrected 'min.counts' mentioned in place of 'min.reads'
# To do:
##  1) Finish script to calculate QC metrics from TEtranscripts info  + multiQC (hisat2) files.
##  2) Allow automated filtering based on either seqmonk QC or TEtranscripts QC (or a mix).
##########################################################################
# Relevant online material:
# 
###########################################################################

# Packages
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(ggplot2,
               meta,
               scuttle,
               cowplot,
               dplyr,
               scater,
               here, 
               stringr)

##################################################################################################################################################################################
# perform_QC_seqmonk_based
##################################################################################################################################################################################

perform_QC_seqmonk_based <- function(input.file = "2_results/colony/HL60/sce.rds",
                                     output.file.dir = '2_results/colony/HL60',
                                     output.fig =  'figures/colony/HL60/QC',
                                     min.counts = 5,
                                     min.cells = 3){


  # Chosen colours for reps (for Al's plots)
  batch_col <- c("rep1" = "#D55E00", "rep2" = "#009E73", "rep3" = "#0072B2")

  #source("src/helpers/library.R")
  suppressWarnings(dir.create(output.fig, recursive = TRUE))
  suppressWarnings(dir.create(dirname(output.file.dir), recursive = TRUE))

  ## ------------------------------------------------------------------------ ##
  cat("\nloading the input file ...\n")
  sce <- readRDS(input.file)


  dupes <- duplicated(rownames(sce))
  sum(dupes)

  # sum(duplicated(scRNA_data$hgnc.TE))
  cat("\nadding QC metrics to sce ...\n")

  # Sub-setting QC metrics
  ## MT
  mito_genes <- which(rowData(sce)$Chromosome %in% 'MT')
  ## TEs
  TEs_sub <- which(rowData(sce)$Type %in% 'TE')
  ## Genes
  all_gene_sub <- which(rowData(sce)$Type %in% 'gene')

  ## Overview of the number of features
  data.frame("Genes" = length(all_gene_sub),
             "TEs" = length(TEs_sub),
             "Mito Genes" = length(mito_genes))



  ##### INITIAL QC! #####################################################################
  ## Subset sce with mito, TEs and Genes for QC

  sce_qc <- scater::addPerCellQC(sce, subsets=list(Mito = mito_genes, 
                                                   TE = TEs_sub, 
                                                   Gene = all_gene_sub))

  ## ------------- mark negative control
  nc <- sce_qc$Treatment == 'Neg'

  ## Extract ColData to make it easier to create custom ggplot figures!

  coldata <- as.data.frame(colData(sce_qc))

  #### Al Plot 1 ########################################################################################
  ## density plots of library sizes across batches
  cat("\nsaving QC plots", output.fig, "...\n")

  #batch_col <- c("rep1" = "#D55E00", "rep2" = "#009E73", "rep3" = "#0072B2")

  ggplot(coldata[!nc,]) + geom_density(aes(total, fill = Batch), alpha = 0.75) +
    scale_fill_manual(values = batch_col) +
    labs(fill = 'Batch', y = 'Density', x = 'Library Size') + theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12,face="bold"))

  ggsave(here(sprintf('%s/Colony_QC1-libsizes.png', output.fig)), width = 7, height = 4)

  #### Al Plot 2 ########################################################################################
  ## plot library size by # of features detected

  coldata$Batch[nc] <- 'Neg'

  ggplot() +
    geom_point(data = coldata, aes(x=detected,
                                   y=total,
                                   size = subsets_Mito_percent,
                                   col = Batch), shape = 1, stroke = 0.4) +
    scale_colour_manual(values =batch_col) +
    labs(x='# of Features Detected',
         y = 'Library Size',
         size = '% Mitochandrial',
        col = 'Group') +
    theme_classic() +
    guides(col = guide_legend(override.aes = list(size = 4)), size = guide_legend(override.aes = list(col = 'grey70'))) +
    expand_limits(x = 12000) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12,face="bold"))

  ggsave(here(sprintf('%s/Colony_QC1-scatter.png', output.fig)), width = 7, height = 4)

  #########################################################################################################################################################
  # Plots (using scater) to produce individual scatter plots for each factor (close to seqmonk output) ###############################################
  head(colData(sce_qc))

  ## Plotted in rows of two
  ### Total counts
  p1 <- plotColData(sce_qc, y = "detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("Number of features detected") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))


  p2 <- plotColData(sce_qc, y = "total", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts (genes + TEs)") + 
    theme(plot.title = element_text(hjust = 0.5))
  ### Mito
  p3 <- plotColData(sce_qc, y = "subsets_Mito_sum", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts in mitochondrial genes") + 
    theme(plot.title = element_text(hjust = 0.5))

  p4 <- plotColData(sce_qc, y = "subsets_Mito_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in mitochondrial genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  ### Genes
  p5 <- plotColData(sce_qc, y = "subsets_Gene_detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("# of genes detected") + 
    theme(plot.title = element_text(hjust = 0.5))

  p6 <- plotColData(sce_qc, y = "subsets_Gene_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in genes") + 
    theme(plot.title = element_text(hjust = 0.5))

  ### TEs
  p7 <- plotColData(sce_qc, y = "subsets_TE_detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("# of TEs detected") + 
    theme(plot.title = element_text(hjust = 0.5))

  p8 <- plotColData(sce_qc, y = "subsets_TE_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in TEs") + 
    theme(plot.title = element_text(hjust = 0.5))


  #### Combine plots and save
  combined1 <- plot_grid(p1, p2,
            p3, p4,
            p5, p6,
            p7, p8,
            ncol =2)

  ggsave(plot = combined1, 
         filename = paste0(output.fig, '/colony_QC1-combined.png'), width = 7, height = 12, bg = 'white')

  ######################################################


  ## Normalised and PCA
  sce_qc_norm <- logNormCounts(sce_qc)
  sce_qc_norm <- runPCA(sce_qc_norm)

  ## Plotting PCA
  plot_grid(
    plotPCA(sce_qc_norm, colour_by = "Batch")+ 
      ggtitle("Rep") + 
      theme(plot.title = element_text(hjust = 0.5)),
  
    plotPCA(sce_qc_norm, colour_by = "Treatment")+ 
      ggtitle(" Treatment") + 
      theme(plot.title = element_text(hjust = 0.5)),
  
    (plotPCA(sce_qc_norm, colour_by = "subsets_Mito_percent") + 
      ggtitle("% of counts in mito genes") + 
      theme(plot.title = element_text(hjust = 0.5)))+ 
      theme(legend.title = element_blank()),

    (plotPCA(sce_qc_norm, colour_by = "subsets_Gene_detected") + 
      ggtitle("# of (unique) genes detected") + 
      theme(plot.title = element_text(hjust = 0.5))) + 
      theme(legend.title = element_blank()),
  ncol = 2)

  ggsave(here(sprintf('%s/colony_QC1-PCA.png', output.fig)), width = 7, height = 7, bg = "white")

  ### Filtering failed samples and low expressed genes
  ############## Secondary QC ###################################################################################################################################
  ## A) Remove samples (columns) failing QC.

  cat("\nRemoving samples that failed QC based on seqmonk results")
  sce_qc2 <- sce_qc[, sce_qc$Passed_seqmonk_QC == "Y"]
  cat("\nRemoved", ncol(sce_qc) - ncol(sce_qc2), "samples that failed Seqmonk. \nLeaving",ncol(sce_qc2),"samples remaining.\n")

  ## B) Remove lowly expressed genes
  #####################
  min_counts <- min.counts
  min_cells <- min.cells
  #####################
  cat("\nFiltering for gene with at least ", min_counts, " reads in ", min_cells, " cells!")

  genes.to.keep <- which(rowSums(counts(sce_qc2) > min_counts) >= min_cells) # Detects which genes/TEs have min_counts expressed in min_cells. Providing a list of IDs and column row numbers to filter and retain.
  sce_qc2 <- sce_qc2[c(genes.to.keep),]

  cat("\nRetained",nrow(sce_qc2), "genes/TEs.")
  cat("\nAfter removing ", nrow(sce_qc)-nrow(sce_qc2)," that didn't meet the above thresholds.\n")
  cat("\nadding QC metrics to sce for samples and gene filtered sce...\n")

  ## Add updated QC info ##################
  # Sub-setting QC metrics
  ## MT
  mito_genes <- which(rowData(sce_qc2)$Chromosome %in% 'MT')
  ## TEs
  TEs_sub <- which(rowData(sce_qc2)$Type %in% 'TE')
  ## Genes
  all_gene_sub <- which(rowData(sce_qc2)$Type %in% 'gene')

  ## Overview of the number of features
  data.frame("Genes" = length(all_gene_sub),
             "TEs" = length(TEs_sub),
             "Mito Genes" = length(mito_genes))
  
  # Remove original QC colData (otherwise it creates duplicates...)
  colData(sce_qc2)[11:ncol(colData(sce_qc2))] <- NULL
  # Add new QC data
  sce_qc2 <- scater::addPerCellQC(sce_qc2, subsets=list(Mito = mito_genes, 
                                                   TE = TEs_sub, 
                                                   Gene = all_gene_sub))

  ## Plotting
  sce_to_plot <- sce_qc2

  ### Total counts
  p1 <- plotColData(sce_to_plot, y = "detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("Number of features detected") + 
    theme(plot.title = element_text(hjust = 0.5))

  p2 <- plotColData(sce_to_plot, y = "total", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts (genes + TEs)") + 
    theme(plot.title = element_text(hjust = 0.5))

  p2.5 <- plotColData(sce_to_plot, y = "total", x = "Batch", colour_by = "Treatment", point_alpha = 0.45) + 
    ggtitle("total counts (genes + TEs)") + 
    theme(plot.title = element_text(hjust = 0.5))
  ### Mito
  p3 <- plotColData(sce_to_plot, y = "subsets_Mito_sum", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts in mitochondrial genes") + 
    theme(plot.title = element_text(hjust = 0.5))

  p4 <- plotColData(sce_to_plot, y = "subsets_Mito_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in mitochondrial genes") + 
    theme(plot.title = element_text(hjust = 0.5))

  p4.5 <- plotColData(sce_to_plot, y = "subsets_Mito_percent", x = "Batch", colour_by = "Treatment", point_alpha = 0.45) + 
    ggtitle("% counts in mitochondrial genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  ### Genes
  p5 <- plotColData(sce_to_plot, y = "subsets_Gene_detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts in genes") + 
    theme(plot.title = element_text(hjust = 0.5))

  p6 <- plotColData(sce_to_plot, y = "subsets_Gene_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in genes") + 
    theme(plot.title = element_text(hjust = 0.5))

  p6.5 <- plotColData(sce_to_plot, y = "subsets_Gene_percent", x = "Batch", colour_by = "Treatment", point_alpha = 0.45) + 
    ggtitle("% counts in genes") + 
    theme(plot.title = element_text(hjust = 0.5))

  ### TEs
  p7 <- plotColData(sce_to_plot, y = "subsets_TE_detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts in TEs") + 
    theme(plot.title = element_text(hjust = 0.5))

  p8 <- plotColData(sce_to_plot, y = "subsets_TE_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in TEs") + 
    theme(plot.title = element_text(hjust = 0.5))

  p8.5 <- plotColData(sce_to_plot, y = "subsets_TE_percent", x = "Batch", colour_by = "Treatment", point_alpha = 0.45) + 
    ggtitle("% counts in TEs") + 
    theme(plot.title = element_text(hjust = 0.5))

  #### Combine plots and save
  plot_grid(p1, p2,p2.5,
            p3, p4,p4.5,
            p5, p6,p6.5,
            p7, p8,p8.5,
            ncol =3)

  ggsave(here(sprintf('%s/colony_QC2-combined.png', output.fig)), width = 9, height = 12, bg = "white")

  ## Normalising
  sce_qc2_norm <- logNormCounts(sce_qc2)
  sce_qc2_norm <- runPCA(sce_qc2_norm)

  ###  Plotting PCA
  plot_grid(
    plotPCA(sce_qc2_norm, colour_by = "Batch")+ 
      ggtitle("Rep") + 
      theme(plot.title = element_text(hjust = 0.5)),
  
    plotPCA(sce_qc2_norm, colour_by = "Treatment")+ 
      ggtitle(" Treatment") + 
      theme(plot.title = element_text(hjust = 0.5)),
  
    (plotPCA(sce_qc2_norm, colour_by = "subsets_Mito_percent") + 
      ggtitle("% of counts in mito genes") + 
      theme(plot.title = element_text(hjust = 0.5)))+ 
      theme(legend.title = element_blank()),
  
    (plotPCA(sce_qc2_norm, colour_by = "subsets_Gene_detected") + 
      ggtitle("# of (unique) genes detected") + 
      theme(plot.title = element_text(hjust = 0.5))) + 
      theme(legend.title = element_blank()),
  ncol = 2)

  ggsave(here(sprintf('%s/colony_QC2-PCA.png', output.fig)), width = 7, height = 7, bg ='white')

  ## Perform PCAs with genes and TEs separately
  sce_qc2_norm <- runPCA(sce_qc2_norm, subset_row = all_gene_sub, name = "PCA_genes")
  sce_qc2_norm <- runPCA(sce_qc2_norm, subset_row = TEs_sub, name = "PCA_TEs")

  plot_grid(
    plotReducedDim(sce_qc2_norm, dimred = "PCA", colour_by = "Treatment")+ 
      ggtitle("PCA 500 variable features (Treatment)") + 
      theme(plot.title = element_text(hjust = 0.5)),
  
    plotReducedDim(sce_qc2_norm, dimred = "PCA", colour_by = "Batch")+ 
      ggtitle("PCA Top 500 variable features (Batch)") + 
      theme(plot.title = element_text(hjust = 0.5)),
  
    plotReducedDim(sce_qc2_norm, dimred = "PCA_genes", colour_by = "Treatment") + 
      ggtitle("PCA genes (Treatment)") + 
      theme(plot.title = element_text(hjust = 0.5)),

    plotReducedDim(sce_qc2_norm, dimred = "PCA_genes", colour_by = "Batch") + 
      ggtitle("PCA genes (Batch)") + 
      theme(plot.title = element_text(hjust = 0.5)),

    plotReducedDim(sce_qc2_norm, dimred = "PCA_TEs", colour_by = "Treatment")+ 
      ggtitle("PCA TEs (Treatment)") + 
      theme(plot.title = element_text(hjust = 0.5)),

    plotReducedDim(sce_qc2_norm, dimred = "PCA_TEs", colour_by = "Batch")+ 
      ggtitle("PCA TEs (Batch)") + 
      theme(plot.title = element_text(hjust = 0.5)),
  ncol = 2)

  ggsave(here(sprintf('%s/colony_QC2-PCA-Genes_TEs_separate.png', output.fig)), width = 7, height = 9, bg = "white")

  ## Compare qc1 and qc2 ########################
  # Obtain variance explained for qc1 and qc2 by
  # 1) Regression
  vars_v1 <- getVarianceExplained(sce_qc_norm, 
                                  variables = c("total", "Batch", "Treatment", "subsets_Gene_detected", "subsets_TE_detected", "subsets_Mito_detected"))
  vars_v2 <- getVarianceExplained(sce_qc2_norm, exprs_values = "logcounts", 
                                  variables = c("total", "Batch", "Treatment", "subsets_Gene_detected", "subsets_TE_detected", "subsets_Mito_detected"))

  ## Plotting
  plot_grid(ncol = 2,
            plotExplanatoryVariables(vars_v1),
            plotExplanatoryVariables(vars_v2))

  ggsave(here(sprintf('%s/colony_QC1.2-explanatoryvariables.png', output.fig)), width = 15, height = 7, bg = "white")


  ################ Save sce files #############################################################################################
  cat("Saving QCd sce data into: ", output.file.dir)

  # Qc1/ pre-filter
  message("Saving pre-filtered sce as: '",output.file.dir,"/sce_qc1.rds'")
  saveRDS(sce_qc, paste0(output.file.dir,'/sce_qc1.rds'))

  # QC2/ post-filter
  message("Saving pre-filtered sce as: '",output.file.dir,"/sce_qc2.rds'")
  saveRDS(sce_qc2, paste0(output.file.dir,'/sce_qc2.rds'))

  message("Minimum number of genes deteced in any sample after QC filtering was:\n", min(sce_qc2$subsets_Gene_detected), " genes")

  if(min(sce_qc2$subsets_Gene_detected) > 10000){
      cat("\n<This is ABOVE the original min TEtranscripts genes measured of 10000")
    } else{
      (cat("\n<This is BELOW the original min TEtranscripts genes measured of 10000") )
  }

}

##################################################################################################################################################################################
# perform_QC_scRNAseq
##################################################################################################################################################################################

perform_QC_scRNAseq <- function(input.file = "./2_results/scNMT/RNAseq/sce.rds",
                                output.file.dir = './2_results/scNMT/RNAseq',
                                output.fig =  'figures/scNMT/RNAseq/QC',
                                min.counts = 5,
                                min.cells = 3,
                                min.cells.percentage = NA, # (0-100) Default applies min.counts in min cells. Any number applied here will then ignore 'min.cells'.
                                thresholds.assay = "counts", # Which assay to apply the thresholds to, default is the original counts. But output will have genes removed across all assays!
                                thresholds.per.batch = FALSE, # If True, the selected thresholds will be applied per batch. Default applies thresholds to all samples, regardless of batch
                                batch.order.col = c("b320" = "#D55E00", "b620" = "#009E73", "b322" = "#0072B2")){
  
  # Chosen colours for reps (for Al's plots)
  #batch_col <- c("b320" = "#D55E00", "b620" = "#009E73", "b322" = "#0072B2")
  batch_col <- batch.order.col
  
  #source("src/helpers/library.R")
  suppressWarnings(dir.create(output.fig, recursive = TRUE))
  suppressWarnings(dir.create(output.file.dir, recursive = TRUE))
  
  ## ------------------------------------------------------------------------ ##
  Line_1 <- paste("\nloading the input file:", input.file,"\n")
  cat(Line_1)
  sce <- readRDS(input.file)
  
  
  dupes <- duplicated(rownames(sce))
  sum(dupes)
  
  # sum(duplicated(scRNA_data$hgnc.TE))
  
  # Sub-setting QC metrics
  ## MT
  mito_genes <- which(rowData(sce)$Chromosome %in% 'MT')
  ## TEs
  TEs_sub <- which(rowData(sce)$Type %in% 'TE')
  ## Genes
  all_gene_sub <- which(rowData(sce)$Type %in% 'gene')
  
  ## Overview of the number of features
  data.frame("Genes" = length(all_gene_sub),
             "TEs" = length(TEs_sub),
             "Mito Genes" = length(mito_genes))
  
  
  
  ##### INITIAL QC! #####################################################################
  
  ## Subset sce with mito, TEs and Genes for QC
  
  sce_qc <- scater::addPerCellQC(sce, subsets=list(Mito = mito_genes, 
                                                   TE = TEs_sub, 
                                                   Gene = all_gene_sub))
  
  
  
  ## ------------- mark negative control
  nc <- sce_qc$Treatment == 'Neg'
  
  ## Extract ColData to make it easier to create custom ggplot figures!
  
  coldata <- as.data.frame(colData(sce_qc))
  
  #### Al Plot 1 ########################################################################################
  ## density plots of library sizes across batches
  Line_2 <- paste("\nsaving QC plots", output.fig, "...\n")
  cat(Line_2)
  
  #batch_col <- c("rep1" = "#D55E00", "rep2" = "#009E73", "rep3" = "#0072B2")
  
  ggplot(coldata[!nc,]) + geom_density(aes(total, fill = Batch), alpha = 0.75) +
    scale_fill_manual(values = batch_col) +
    labs(fill = 'Batch', y = 'Density', x = 'Library Size') + theme_bw() +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12,face="bold"))
  
  ggsave(here(sprintf('%s/Colony_QC1-libsizes.png', output.fig)), width = 7, height = 4)
  
  
  #### Al Plot 2 ########################################################################################
  ## plot library size by # of features detected
  
  coldata$Batch[nc] <- 'Neg'
  
  ggplot() +
    geom_point(data = coldata, aes(x=detected,
                                   y=total,
                                   size = subsets_Mito_percent,
                                   col = Batch), shape = 1, stroke = 0.4) +
    scale_colour_manual(values =batch_col) +
    labs(x='# of Features Detected',
         y = 'Library Size',
         size = '% Mitochandrial',
         col = 'Group') +
    theme_classic() +
    guides(col = guide_legend(override.aes = list(size = 4)), size = guide_legend(override.aes = list(col = 'grey70'))) +
    expand_limits(x = 12000) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12,face="bold"))
  
  ggsave(here(sprintf('%s/Colony_QC1-scatter.png', output.fig)), width = 7, height = 4)
  
  
  
  
  #########################################################################################################################################################
  # Plots (using scater) to produce individual scatter plots for each factor (close to seqmonk output) ###############################################
  
  head(colData(sce_qc))
  colData(sce_qc)$Batch <- factor(colData(sce_qc)$Batch, levels = names(batch.order.col))
  
  ## Plotted in rows of two
  ### Total counts
  p1 <- plotColData(sce_qc, y = "detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("Number of features detected") + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  p2 <- plotColData(sce_qc, y = "total", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts (genes + TEs)") + 
    theme(plot.title = element_text(hjust = 0.5))
  ### Mito
  p3 <- plotColData(sce_qc, y = "subsets_Mito_sum", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts in mitochondrial genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p4 <- plotColData(sce_qc, y = "subsets_Mito_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in mitochondrial genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  ### Genes
  p5 <- plotColData(sce_qc, y = "subsets_Gene_detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("# of genes detected") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p6 <- plotColData(sce_qc, y = "subsets_Gene_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  ### TEs
  p7 <- plotColData(sce_qc, y = "subsets_TE_detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("# of TEs detected") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p8 <- plotColData(sce_qc, y = "subsets_TE_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in TEs") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  #### Combine plots and save
  
  combined1 <- plot_grid(p1, p2,
                         p3, p4,
                         p5, p6,
                         p7, p8,
                         ncol =2)
  
  
  ggsave(plot = combined1, 
         filename = paste0(output.fig, '/colony_QC1-combined.png'), width = 7, height = 12, bg = 'white')
  
  ######################################################
  ## Normalised and PCA
  sce_qc_norm <- logNormCounts(sce_qc)
  sce_qc_norm <- runPCA(sce_qc_norm)
  
  ## Plotting PCA
  plot_grid(
    plotPCA(sce_qc_norm, colour_by = "Batch") + 
      ggtitle("Rep") + 
      theme(plot.title = element_text(hjust = 0.5)),
    
    plotPCA(sce_qc_norm, colour_by = "Treatment") + 
      ggtitle(" Treatment") + 
      theme(plot.title = element_text(hjust = 0.5)),
    
    (plotPCA(sce_qc_norm, colour_by = "subsets_Mito_percent") + 
      ggtitle("% of counts in mito genes") + 
      theme(plot.title = element_text(hjust = 0.5))) + 
      theme(legend.title = element_blank()),
    
    (plotPCA(sce_qc_norm, colour_by = "subsets_Gene_detected") + 
      ggtitle("# of (unique) genes detected") + 
      theme(plot.title = element_text(hjust = 0.5))) + 
      theme(legend.title = element_blank()),
    ncol = 2)
  
  ggsave(here(sprintf('%s/colony_QC1-PCA.png', output.fig)), width = 7, height = 7, bg = "white")
  
  ### Filtering failed samples and low expressed genes
  ############## Secondary QC ###################################################################################################################################
  ## A) Remove samples (columns) failing QC.
  Line_3 <- paste("\nRemoving samples that failed QC based on seqmonk results")
  cat(Line_3)
  
  sce_qc2 <- sce_qc[, sce_qc$Pass_scNMT_QC == "Y"|TRUE]
  
  Line_4 <- paste("\nRemoved", ncol(sce_qc) - ncol(sce_qc2), "samples that failed Seqmonk. \nLeaving",ncol(sce_qc2),"samples remaining.\n")
  cat(Line_4)
  
  ## B) Remove lowly expressed genes
  ### by NUMBER of 'min.cells'. ###############
  
  if(is.na(min.cells.percentage) == TRUE){
    
    ### By min.counts in min.cells
    #####################
    min_counts <- min.counts
    min_cells <- min.cells
    #####################
    ### Across all cells ###############
    if(thresholds.per.batch == FALSE){
      cat("Applying filtering thresholds across all cells (regardless of batch)")
      
      Line_5 <- paste("\nFiltering for gene with at least ", min_counts, " reads in ", min_cells, " cells!")
      cat(Line_5)
      #####################
      genes.to.keep <- which(rowSums(assay(sce_qc2, thresholds.assay) > min_counts) >= min_cells) # Detects which genes/TEs have min_counts expressed in min_cells. Providing a list of IDs and column row numbers to filter and retain.
      sce_qc2 <- sce_qc2[c(genes.to.keep),]
      #####################
    } else if(thresholds.per.batch == TRUE){
      
    ### Applied per Batch ###############
      batch_thresholds <- 
        table(sce_qc2$Batch) %>% 
        as.data.frame(.) %>% 
        rename("batch" = Var1, "n_samples" = Freq) %>% 
        mutate("min_counts" = min_counts,
               "min_cells" = ceiling(min_cells),
               batch = as.character(batch))
      
      Line_5 <- paste("\nRequiring the following minimal counts and cells per batch (see 'Passing_QC_thresholds_per_batch.csv' for stored info) ")
      cat(Line_5)
      batch_thresholds
      
      genes.to.keep.list <- list()
      n_features_per_batch <- c()
      for(batch in batch_thresholds$batch) {
        
        assay_to_check <- assay(sce_qc2[,sce_qc2$Batch == batch], thresholds.assay) # Filter for thresholds for specific batch
        batch_thresholds_subset <- filter(batch_thresholds, batch == batch)
        
        genes.to.keep.list[[batch]] <- names(which(rowSums(assay_to_check > min_counts) >= batch_thresholds_subset$min_cells))
        n_features_per_batch <- c(n_features_per_batch, length(names(which(rowSums(assay_to_check > min_counts) >= batch_thresholds_subset$min_cells))))
      }
      
      genes_in_common <- Reduce(intersect, genes.to.keep.list)
      
      batch_thresholds$Features_passing_threshold <- n_features_per_batch
      batch_thresholds$Features_overlap_retained <- length(genes_in_common)
      
      sce_qc2 <- sce_qc2[c(genes_in_common),]
      
      cat("\nSaving filtered info: ",paste0(output.file.dir,"/Passing_QC_thresholds_per_batch.csv"))
      write.csv(batch_thresholds, file = paste0(output.file.dir,"/Passing_QC_thresholds_per_batch.csv"), row.names = F)
      batch_thresholds
    }
    
    
  ### by PERCENTAGE of min cells ############
    
  }  else if(is.na(min.cells.percentage) == FALSE){
    
    ### By min.counts in % of cells - rounded up to the nearest integer.
    #####################
    min_counts <- min.counts
    min_cells_percentage <- min.cells.percentage
    ########################
    
    ### Across all cells ###############
    if(thresholds.per.batch == FALSE){
      cat("Applying filtering thresholds across all cells (regardless of batch)")
      
      min_cells <- ceiling(ncol(sce_qc2)*(min_cells_percentage/100))
      Line_5 <- paste("\nFiltering for gene with at least ", min_counts, " reads in ", min_cells_percentage, "% of cells:\n",min_cells, " out of ", ncol(sce_qc2), " cells")
      cat(Line_5)
      #####################
      genes.to.keep <- which(rowSums(assay(sce_qc2, thresholds.assay) > min_counts) >= min_cells) # Detects which genes/TEs have min_counts expressed in min_cells. Providing a list of IDs and column row numbers to filter and retain.
      sce_qc2 <- sce_qc2[c(genes.to.keep),]
      
    ### Applied per Batch ###############
    } else if(thresholds.per.batch == TRUE){
      cat("\nApplying filtering thresholds across cells per batch")
      
      batch_thresholds <- 
        table(sce_qc2$Batch) %>% 
        as.data.frame(.) %>% 
        rename("batch" = Var1, "n_samples" = Freq) %>% 
        mutate("perc_filter" = min_cells_percentage/100,
               "min_cells" = ceiling(n_samples*perc_filter),
               batch = as.character(batch))
      
      Line_5 <- paste("\nRequiring the following minimal counts and cells per batch (see 'Passing_QC_thresholds_per_batch.csv' for stored info) ")
      cat(Line_5)
      batch_thresholds
      
      genes.to.keep.list <- list()
      n_features_per_batch <- c()
      for(batch in batch_thresholds$batch) {
        
        assay_to_check <- assay(sce_qc2[,sce_qc2$Batch == batch], thresholds.assay) # Filter for thresholds for specific batch
        batch_thresholds_subset <- filter(batch_thresholds, batch == batch)
        
        genes.to.keep.list[[batch]] <- names(which(rowSums(assay_to_check > min_counts) >= batch_thresholds_subset$min_cells))
        n_features_per_batch <- c(n_features_per_batch, length(names(which(rowSums(assay_to_check > min_counts) >= batch_thresholds_subset$min_cells))))
        
      }
      
      genes_in_common <- Reduce(intersect, genes.to.keep.list)
      
      batch_thresholds$Features_passing_threshold <- n_features_per_batch
      batch_thresholds$Features_overlap_retained <- length(genes_in_common)
      
      sce_qc2 <- sce_qc2[c(genes_in_common),]
      
      cat("\nSaving filtered info: ",paste0(output.file.dir,"/Passing_QC_thresholds_per_batch.csv"))
      write.csv(batch_thresholds, file = paste0(output.file.dir,"/Passing_QC_thresholds_per_batch.csv"), row.names = F)
      batch_thresholds
    }
  }
  
  # Counting genes vs TEs retained
  qc1_genes <- sum(str_count(rownames(sce_qc), "ENSG"))
  qc1_TEs <- sum(!str_count(rownames(sce_qc), "ENSG"))
  qc2_genes <- sum(str_count(rownames(sce_qc2), "ENSG"))
  qc2_TEs <- sum(!str_count(rownames(sce_qc2), "ENSG"))
  # Report above numbers
  Line_6 <- paste("\nRetained ", (qc2_genes+qc2_TEs), " features \n(",qc2_genes," genes / ",qc2_TEs," TEs)", sep = "")
  Line_7 <- paste("\nAfter removing ", nrow(sce_qc)-nrow(sce_qc2)," that didn't meet the above thresholds.\n",
                  "(",qc1_genes - qc2_genes, " genes removed / ", qc1_TEs - qc2_TEs , " TEs removed)", sep ="")
  Line_8 <- paste("\nadding QC metrics to sce for samples and gene filtered sce...\n")
  cat(Line_6)
  cat(Line_7)
  cat(Line_8)

  ## Add updated QC info ##################
  
  # Sub-setting QC metrics
  ## MT
  mito_genes <- which(rowData(sce_qc2)$Chromosome %in% 'MT')
  ## TEs
  TEs_sub <- which(rowData(sce_qc2)$Type %in% 'TE')
  ## Genes
  all_gene_sub <- which(rowData(sce_qc2)$Type %in% 'gene')
  
  ## Overview of the number of features
  data.frame("Genes" = length(all_gene_sub),
             "TEs" = length(TEs_sub),
             "Mito Genes" = length(mito_genes))
  
  
  # Remove original QC colData (otherwise it creates duplicates...)
  colData(sce_qc2)[11:ncol(colData(sce_qc2))] <- NULL
  # Add new QC data
  sce_qc2 <- scater::addPerCellQC(sce_qc2, subsets=list(Mito = mito_genes, 
                                                        TE = TEs_sub, 
                                                        Gene = all_gene_sub))
  
  ## Plotting
  sce_to_plot <- sce_qc2
  
  ### Total counts
  p1 <- plotColData(sce_to_plot, y = "detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("Number of features detected") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  
  p2 <- plotColData(sce_to_plot, y = "total", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts (genes + TEs)") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p2.5 <- plotColData(sce_to_plot, y = "total", x = "Batch", colour_by = "Treatment", point_alpha = 0.45) + 
    ggtitle("total counts (genes + TEs)") + 
    theme(plot.title = element_text(hjust = 0.5))
  ### Mito
  p3 <- plotColData(sce_to_plot, y = "subsets_Mito_sum", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts in mitochondrial genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p4 <- plotColData(sce_to_plot, y = "subsets_Mito_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in mitochondrial genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p4.5 <- plotColData(sce_to_plot, y = "subsets_Mito_percent", x = "Batch", colour_by = "Treatment", point_alpha = 0.45) + 
    ggtitle("% counts in mitochondrial genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  ### Genes
  p5 <- plotColData(sce_to_plot, y = "subsets_Gene_detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts in genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p6 <- plotColData(sce_to_plot, y = "subsets_Gene_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p6.5 <- plotColData(sce_to_plot, y = "subsets_Gene_percent", x = "Batch", colour_by = "Treatment", point_alpha = 0.45) + 
    ggtitle("% counts in genes") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  ### TEs
  p7 <- plotColData(sce_to_plot, y = "subsets_TE_detected", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("total counts in TEs") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p8 <- plotColData(sce_to_plot, y = "subsets_TE_percent", x = "Batch", colour_by = "Batch", point_alpha = 0.45) + 
    ggtitle("% counts in TEs") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  p8.5 <- plotColData(sce_to_plot, y = "subsets_TE_percent", x = "Batch", colour_by = "Treatment", point_alpha = 0.45) + 
    ggtitle("% counts in TEs") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  
  #### Combine plots and save
  
  plot_grid(p1, p2,p2.5,
            p3, p4,p4.5,
            p5, p6,p6.5,
            p7, p8,p8.5,
            ncol =3)
  
  ggsave(here(sprintf('%s/colony_QC2-combined.png', output.fig)), width = 9, height = 12, bg = "white")
  
  ## Normalising
  sce_qc2_norm <- logNormCounts(sce_qc2)
  sce_qc2_norm <- runPCA(sce_qc2_norm)
  
  ###  Plotting PCA
  plot_grid(
    
    plotPCA(sce_qc2_norm, colour_by = "Batch")+ 
      ggtitle("Rep") + 
      theme(plot.title = element_text(hjust = 0.5)),
    
    plotPCA(sce_qc2_norm, colour_by = "Treatment")+ 
      ggtitle(" Treatment") + 
      theme(plot.title = element_text(hjust = 0.5)),
    
    (plotPCA(sce_qc2_norm, colour_by = "subsets_Mito_percent") + 
      ggtitle("% of counts in mito genes") + 
      theme(plot.title = element_text(hjust = 0.5)))+ 
      theme(legend.title = element_blank()),
    
    (plotPCA(sce_qc2_norm, colour_by = "subsets_Gene_detected") + 
      ggtitle("# of (unique) genes detected") + 
      theme(plot.title = element_text(hjust = 0.5))) + 
      theme(legend.title = element_blank()),
    ncol = 2)
  
  ggsave(here(sprintf('%s/colony_QC2-PCA.png', output.fig)), width = 7, height = 7, bg ='white')
  
  ## Perform PCAs with genes and TEs separately
  # all_gene_sub
  
  sce_qc2_norm <- runPCA(sce_qc2_norm, subset_row = all_gene_sub, name = "PCA_genes")
  sce_qc2_norm <- runPCA(sce_qc2_norm, subset_row = TEs_sub, name = "PCA_TEs")
  
  plot_grid(
    plotReducedDim(sce_qc2_norm, dimred = "PCA", colour_by = "Treatment")+ 
      ggtitle("PCA 500 variable features (Treatment)") + 
      theme(plot.title = element_text(hjust = 0.5)),
    
    plotReducedDim(sce_qc2_norm, dimred = "PCA", colour_by = "Batch")+ 
      ggtitle("PCA Top 500 variable features (Batch)") + 
      theme(plot.title = element_text(hjust = 0.5)),
    
    plotReducedDim(sce_qc2_norm, dimred = "PCA_genes", colour_by = "Treatment") + 
      ggtitle("PCA genes (Treatment)") + 
      theme(plot.title = element_text(hjust = 0.5)),
    
    plotReducedDim(sce_qc2_norm, dimred = "PCA_genes", colour_by = "Batch") + 
      ggtitle("PCA genes (Batch)") + 
      theme(plot.title = element_text(hjust = 0.5)),
    
    plotReducedDim(sce_qc2_norm, dimred = "PCA_TEs", colour_by = "Treatment")+ 
      ggtitle("PCA TEs (Treatment)") + 
      theme(plot.title = element_text(hjust = 0.5)),
    
    plotReducedDim(sce_qc2_norm, dimred = "PCA_TEs", colour_by = "Batch")+ 
      ggtitle("PCA TEs (Batch)") + 
      theme(plot.title = element_text(hjust = 0.5)),
    ncol = 2)
  
  ggsave(here(sprintf('%s/colony_QC2-PCA-Genes_TEs_separate.png', output.fig)), width = 7, height = 9, bg = "white")
  
  ## Compare qc1 and qc2 ########################
  # Obtain variance explained for qc1 and qc2 by
  # 1) Regression
  vars_v1 <- getVarianceExplained(sce_qc_norm, 
                                  variables = c("total", "Batch", "Treatment", "subsets_Gene_detected", "subsets_TE_detected", "subsets_Mito_detected"))
  vars_v2 <- getVarianceExplained(sce_qc2_norm, exprs_values = "logcounts", 
                                  variables = c("total", "Batch", "Treatment", "subsets_Gene_detected", "subsets_TE_detected", "subsets_Mito_detected"))
  
  ## Plotting
  plot_grid(ncol = 2,
            plotExplanatoryVariables(vars_v1),
            plotExplanatoryVariables(vars_v2))
  
  ggsave(here(sprintf('%s/colony_QC1.2-explanatoryvariables.png', output.fig)), width = 15, height = 7, bg = "white")
  
  ################ Save sce files #############################################################################################
  Line_9 <- paste("Saving QCd sce data into: ", output.file.dir)
  cat(Line_9)
  
  
  # Qc1/ pre-filter
  #output.file = '2_results/colony/HL60_updated/sce_qc.rds'
  message("Saving pre-filtered sce as: '",output.file.dir,"/sce_qc1.rds'")
  saveRDS(sce_qc, paste0(output.file.dir,'/sce_qc1.rds'))
  
  # QC2/ post-filter
  #output.file = '2_results/colony/HL60_updated/sce_qc2.rds'
  message("Saving pre-filtered sce as: '",output.file.dir,"/sce_qc2.rds'")
  saveRDS(sce_qc2, paste0(output.file.dir,'/sce_qc2.rds'))
  
  Line_10 <- paste("\nMinimum number of genes deteced in any sample after QC filtering was:\n", min(sce_qc2$subsets_Gene_detected), " genes")
  cat(Line_10)
  
  # Save relevant filtering info into text file:
  writeLines(c(Line_1,Line_2,Line_3,Line_4,Line_5,Line_6,Line_7,Line_8,Line_9,Line_10), 
             con = paste0(output.file.dir,"/QC_info.txt"), sep = "")
  
  if(min(sce_qc2$subsets_Gene_detected) > 10000){
    cat("\n<This is ABOVE the original min TEtranscripts genes measured of 10000")
  } else {
    cat("\n<This is BELOW the original min TEtranscripts genes measured of 10000")
  }
  
}





