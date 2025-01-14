##########################################################################
# Purpose: R Script to perform per batch normalisation with scran/scater package
# Output: Produces four text files of samples passing QC:
## i) list of sample IDs passing QC in a single column
#
# Date: 25.August.22
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: This was originally developed as part of a package build running on unix/HPC by Al, which would allow shell style submission.
# To do: Automate running of visualisation script (add function at end).
## option for perBatch and all batches to normalise by:
## scuttle:logNormCounts
## scran - preCluster and pooled
## edgeR - RLE
## Make option for scuttle not scale by gene
## Add plot of norm factor calculated/used vs total counts (library size)

##########################################################################

# Load packages
library(ggplot2)
#library(meta)
#library(scuttle)
library(cowplot)
library(dplyr)
library(scater)
library(scran)
library(edgeR)


# Test options
#input.file = "./2_results/colony/HL60/sce_qc2.rds"
#output.file = "./2_results/colony/HL60/Batch_Correction/sce_qc2_PerBatchNormPCA.rds"
#force = FALSE

########################################################################################################
## RLE
########################################################################################################

per_batch_normalisation_edge_RLE <- function(input.file = "./2_results/scNMT/RNAseq/sce_qc2.rds",
                                    output.file = "./2_results/scNMT/RNAseq/Batch_Correction/sce_qc2_PerBatchNormPCA-RLE.rds",
                                    output.fig.dir = "./figures/scNMT/RNAseq/Batch_Correction/Normalisation/",
                                    force = FALSE,
                                    scale.genes = TRUE,
                                    per.batch = TRUE){ # Whether to apply RLE to all samples or split them into individual batches first then apply RLE and runPCA then combine back after.
  
  
  # TO DELETE
  #input.file = "./2_results/scNMT/RNAseq/sce_qc2.rds"
  #output.file = "./2_results/scNMT/RNAseq/Batch_Correction/sce_qc2_PerBatchNormPCA-RLE.rds"
  #output.fig.dir = "./figures/scNMT/RNAseq/Batch_Correction/Normalisation/"
  #  force = FALSE
  #  scale.genes = TRUE
  #  per.batch = TRUE
  
  
  # Check existence and create dir if necessary
  
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  
  
  # Read in file
  cat("\nloading the input file:", input.file)
  sce_qc <- readRDS(input.file)
  unique(sce_qc$Batch)
  
  ## ------------- Multi-Batch Normalisation
  cat("\nperforming per-batch normalisation ...\n")
  cat("\nwith post normalisation PCA and umap calculations\n")
  sce.list <- list()
  
  
  
  if(per.batch == FALSE){
    # Calculate RLE with edgeR
    dgeSeq = edgeR::DGEList(counts=assay(sce_qc, 'counts'))
    dgeSeq2 = edgeR::calcNormFactors(dgeSeq, method="RLE")
    # Save Norm/Size.Factor
    # sce_qc$sizeFactor.RLE <- dgeSeq2$samples$norm.factors
    sce_qc$sizeFactor <- as.vector(dgeSeq2$samples$norm.factors)
    sce_qc$sizeFactor.RLE <- as.vector(dgeSeq2$samples$norm.factors)
    
    # LogNormCounts using RLE
    sce_qc2 <- logNormCounts(sce_qc) #, center.size.factors = F) # This automatically uses the sizeFactor column and stores it in 'logcounts'
    
    # Log2 CPM using normalised factors and just using 'library size' via scuttle::calculateCPM() vs RLE + edgeR::cpm()
    assay(sce_qc2, "log.CPM.RLE") = log2(calculateCPM(sce_qc2, size.factors = as.vector(sce_qc2$sizeFactor)) + 1)
    assay(sce_qc2, "log.CPM.libsize") = log2(calculateCPM(sce_qc2, size.factors = as.vector(sce_qc2$total)) + 1)
    assay(sce_qc2, "log.CPM.edgeR.RLE") = edgeR::cpm(dgeSeq2, prior.count = 2, log = T)
    
    
    # Save scatter plots of Norm Factors vs library size.
    p1 <- as.data.frame(colData(sce_qc2)) %>% ggplot(aes(x=total, y =sizeFactor.RLE, color = Batch, shape = Treatment)) + geom_point() + coord_trans(x="log2") + theme_light() + xlab("Library size")
    p2 <- as.data.frame(colData(sce_qc2)) %>% ggplot(aes(x=total, y =sizeFactor.RLE, color = Treatment, shape = Batch)) + geom_point() + coord_trans(x="log2") + theme_light() + xlab("Library size")+
      scale_colour_manual(values = c("Unt" = "green", "AZA" = "blue", "DAC" = "red"))
    plot_grid(p1, p2, labels = c('A', 'B'))
    ggsave(filename = paste0(output.fig.dir, "SizeFactor.RLE_vs_lib.size_perbatchFalse.png"), width = 18, height = 6)
    
    # runPCA on individual batches and merge into a single dataframe
    
    Batches <- unique(sce_qc2$Batch)
    for (Batch in unique(sce_qc2$Batch)){
      set.seed(405)
      sce_i <- sce_qc2[,sce_qc2$Batch == Batch]
      sce_i <- runPCA(sce_i, ncomponents = 4, name = 'PCA.batch')
      sce_i <- runUMAP(sce_i, ncomponents = 4, name = 'UMAP.batch')
      sce.list[[Batch]] <- sce_i
    }
    
    sce.batch <- Reduce('cbind', sce.list)  
    
  } else if(per.batch == TRUE){
    
    Batches <- unique(sce_qc$Batch)
    
    for (Batch in unique(sce_qc$Batch))
    {
      set.seed(405)
      
      sce_i <- sce_qc[,sce_qc$Batch == Batch]
      
      
      dgeSeq = edgeR::DGEList(counts=assay(sce_i, 'counts'))
      dgeSeq = edgeR::calcNormFactors(dgeSeq, method="RLE")
      
      sce_i$sizeFactor <- as.vector(dgeSeq$samples$norm.factors)
      sce_i$sizeFactor.RLE <- as.vector(dgeSeq$samples$norm.factors) # Saved in two slots to make it easier comparison for confirming the correct size.factor was applied.
      
      sce_i <- logNormCounts(sce_i)
      
      # Log2 CPM using normalised factors and just using 'library size' via scuttle::calculateCPM() vs RLE + edgeR::cpm()
      assay(sce_i, "log.CPM.RLE") = log2(calculateCPM(sce_i, size.factors = as.vector(sce_i$sizeFactor)) + 1)
      assay(sce_i, "log.CPM.libsize") = log2(calculateCPM(sce_i, size.factors = as.vector(sce_i$total)) + 1)
      assay(sce_i, "log.CPM.edgeR.RLE") = edgeR::cpm(dgeSeq, prior.count = 2, log = T)
      
      sce_i <- runPCA(sce_i, ncomponents = 4, name = 'PCA.batch')
      sce_i <- runUMAP(sce_i, ncomponents = 4, name = 'UMAP.batch')
      sce.list[[Batch]] <- sce_i
    }
    sce.batch <- Reduce('cbind', sce.list)
    
    p1 <- as.data.frame(colData(sce.batch)) %>% ggplot(aes(x=total, y =sizeFactor.RLE, color = Batch, shape = Treatment)) + geom_point() + coord_trans(x="log2") + theme_light() + xlab("Library size")
    p2 <- as.data.frame(colData(sce.batch)) %>% ggplot(aes(x=total, y =sizeFactor.RLE, color = Treatment, shape = Batch)) + geom_point() + coord_trans(x="log2") + theme_light() + xlab("Library size")
    plot_grid(p1, p2, labels = c('A', 'B'))
    ggsave(filename = paste0(output.fig.dir, "SizeFactor.RLE_vs_lib.size_perbatchTRUE.png"), width = 18, height = 6)
  }
  
  # runPCA with batches together post normalisation
  sce.batch <- runPCA(sce.batch, ncomponents = 4, name = 'PCA.batch.united', exprs_values = "logcounts")
  
  #
  cat("writing output file ...\n")
  message(output.file)
  saveRDS(sce.batch, output.file)
  
  cat("\nFinished normalisation script. Now run visualisaiton script.\n")
  
}

########################################################################################################
## Pooled_Deconvolution
########################################################################################################

per_batch_normalisation_Pooled_Deconvolution <- function(input.file = "./2_results/scNMT/RNAseq/sce_qc2.rds",
                                                         output.file = "./2_results/scNMT/RNAseq/Batch_Correction_scran/sce_qc2_Norm-PoolDeconv.rds",
                                                         output.fig.dir = "./figures/scNMT/RNAseq/Batch_Correction/Normalisation/",
                                                         apply.quickCluster = FALSE, # This is probably not too useful with this dataset as many clusters simply identify groups within each batch. But maybe that would at least aid preservation of similar cells within those groups!?
                                                         quickCluster.method = "igraph", # Default is igraph. Alternative is 'hclust' or 'NONE'. If none. min.size will be ignored (no need to change it).
                                                         quickCluster.min.size = 20,
                                                         Pooled.window.sizes = seq(21,101,5),
                                                         perBatch = FALSE, # But runPCA is automatically still performed for batches individually
                                                         force = FALSE){
  
  
  ## To delete
  
#  input.file = "./2_results/scNMT/RNAseq/sce_qc2.rds"
#  output.file = "./2_results/scNMT/RNAseq/Batch_Correction_scran/Normalisation/sce_qc2_Norm-PoolDeconv.rds"
#  apply.quickCluster = FALSE # This is probably not too useful with this dataset as many clusters simply identify groups within each batch. But maybe that would at least aid preservation of similar cells within those groups!?
#  quickCluster.method = "igraph" # Default is igraph. Alternative is 'hclust' or 'NONE'. If none. min.size will be ignored (no need to change it).
#  quickCluster.min.size = 20
#  Pooled.window.sizes = seq(21,101,5)
#  perBatch = FALSE # But runPCA is automatically still performed for batches individually
#  force = FALSE
  
  # This applies the pooled-based deconvolution approach from Lun et al 2016. Which was first written in the scran pacakge and migrated to the scuttle package.
  # As this relies on pooling groups of cells (default was pools of 100 cells) this is not written to allow splitting of batches for individual normalisation.
  # The PCA after normalisation are calculated on the batches independently to identify the most variable genes/TEs specific to their own batches.
  # A separate PCA is run on newly 'normalised' data to find the most variable genes - most likely this will pick up genes that explain differences in batches. Hence the requirement of additional batch correction methods...
  
  
  # Check existence and create dir if necessary
  
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  
  
  # Read in file
  cat("\nloading the input file:", input.file)
  sce_qc <- readRDS(input.file)
  unique(sce_qc$Batch)
  
  ## ------------- Multi-Batch Normalisation
  cat("\nperforming per-batch normalisation ...\n")
  cat("\nwith post normalisation PCA and umap calculations\n")
  sce.list <- list()
  
  
  # Computing the size factors
  # Using pooled with pre-clustering and scuttle/scran
  ## Perform Clustering of cells prior for pooled devonvolution approach.
  library(scran)
  if(apply.quickCluster == TRUE){
    
  preclusters <- quickCluster(sce_qc, min.size = quickCluster.min.size, method = quickCluster.method) # Modified to 20 from default 50 to allow for 1/3 of a batch (roughly one treatment group) to be the minimal cluster size. Arguably this could slightly smaller...
  table(preclusters)
  sce_qc_scuttle <- computePooledFactors(sce_qc, sizes = Pooled.window.sizes, clusters=preclusters)
  
  sce_qc_scuttle$sizeFactor.deconv <- sce_qc_scuttle$sizeFactor # Make a copy for easier confirmation of which factor was applied later on as default options tend to use 'sizeFactor' slot.
  sce_qc_scuttle$sizeFactor.deconv2 <- scuttle::pooledSizeFactors(sce_qc_scuttle, scaling=sce_qc_scuttle$sizeFactor) # Applying it a second time with first round scaling supplied can aid accuracy.
  
  sce_qc_scuttle <- logNormCounts(sce_qc_scuttle, size.factors =sce_qc_scuttle$sizeFactor.deconv)
  sce_qc_scuttle <- logNormCounts(sce_qc_scuttle, size.factors =sce_qc_scuttle$sizeFactor.deconv2, name = "logcounts.sf.deconv2")
  
  } else {
  
  sce_qc_scuttle <- computePooledFactors(sce_qc, sizes = Pooled.window.sizes)
  sce_qc_scuttle$sizeFactor.deconv <- sce_qc_scuttle$sizeFactor 
  sce_qc_scuttle$sizeFactor.deconv2 <- scuttle::pooledSizeFactors(sce_qc_scuttle, scaling=sce_qc_scuttle$sizeFactor) # Applying it a second time with first round scaling supplied can aid accuracy.
  
  sce_qc_scuttle <- logNormCounts(sce_qc_scuttle, size.factors =sce_qc_scuttle$sizeFactor.deconv)
  sce_qc_scuttle <- logNormCounts(sce_qc_scuttle, size.factors =sce_qc_scuttle$sizeFactor.deconv2, name = "logcounts.sf.deconv2")
  
  }
  
  
  p1 <- as.data.frame(colData(sce_qc_scuttle)) %>% ggplot(aes(x=total, y =sizeFactor.deconv, color = Batch, shape = Treatment)) + geom_point() + coord_trans(x="log2", y = "log2") + theme_light() + xlab("Library size") 
  p2 <- as.data.frame(colData(sce_qc_scuttle)) %>% ggplot(aes(x=total, y =sizeFactor.deconv, color = Treatment, shape = Batch)) + geom_point() + coord_trans(x="log2", y = "log2") + theme_light() + xlab("Library size")+
    scale_colour_manual(values = c("Unt" = "green", "AZA" = "blue", "DAC" = "red"))
  plot_grid(p1, p2, labels = c('A', 'B'))
  ggsave(filename = paste0(output.fig.dir, "SizeFactor.pool.deconv_vs_lib.sizeE.png"), width = 18, height = 6)
  
  p1 <- as.data.frame(colData(sce_qc_scuttle)) %>% ggplot(aes(x=total, y =sizeFactor.deconv2, color = Batch, shape = Treatment)) + geom_point() + coord_trans(x="log2", y = "log2") + theme_light() + xlab("Library size")
  p2 <- as.data.frame(colData(sce_qc_scuttle)) %>% ggplot(aes(x=total, y =sizeFactor.deconv2, color = Treatment, shape = Batch)) + geom_point() + coord_trans(x="log2", y = "log2") + theme_light() + xlab("Library size")+
    scale_colour_manual(values = c("Unt" = "green", "AZA" = "blue", "DAC" = "red"))
  plot_grid(p1, p2, labels = c('A', 'B'))
  ggsave(filename = paste0(output.fig.dir, "SizeFactor.pool.deconv2_vs_lib.size.png"), width = 18, height = 6)
  
  
  # runPCA on individual batches and merge into a single dataframe
  
  Batches <- unique(sce_qc_scuttle$Batch)
  for (Batch in unique(sce_qc_scuttle$Batch)){
    set.seed(405)
    sce_i <- sce_qc_scuttle[,sce_qc_scuttle$Batch == Batch]
    Genes_Index <- which(rowData(sce_i)$Type %in% 'gene')
    TEs_Index <- which(rowData(sce_i)$Type %in% 'TE')
    
    sce_i <- runPCA(sce_i, ncomponents = 4, name = 'PCA.batch', exprs_values = "logcounts")
    sce_i <- runUMAP(sce_i, ncomponents = 4, name = 'UMAP.batch-Genes', exprs_values = "logcounts", subset_row = Genes_Index)
    sce_i <- runUMAP(sce_i, ncomponents = 4, name = 'UMAP.batch-TEs', exprs_values = "logcounts", subset_row = TEs_Index)
    
    sce_i <- runPCA(sce_i, ncomponents = 4, name = 'PCA.batch.deconv2', exprs_values = "logcounts.sf.deconv2")
    sce_i <- runUMAP(sce_i, ncomponents = 4, name = 'UMAP.batch.deconv2-Genes', exprs_values = "logcounts.sf.deconv2", subset_row = Genes_Index)
    sce_i <- runUMAP(sce_i, ncomponents = 4, name = 'UMAP.batch.deconv2-TEs', exprs_values = "logcounts.sf.deconv2", subset_row = TEs_Index)
    
    sce.list[[Batch]] <- sce_i
  }
  
  sce.batch <- Reduce('cbind', sce.list)  
  
  # runPCA with batches together post normalisation
  sce.batch <- runPCA(sce.batch, ncomponents = 4, name = 'PCA.batch.united', exprs_values = "logcounts")
  sce.batch <- runPCA(sce.batch, ncomponents = 4, name = 'PCA.batch.united.sf.deconv2', exprs_values = "logcounts.sf.deconv2")

  # TODO variance explained by each PC
  cat("writing output file ...\n")
  message(output.file)
  saveRDS(sce.batch, output.file)
  
  cat("\nFinished normalisation script. Now run visualisaiton script.\n")
  
}





########################################################################################################
## per_batch_normalisation
########################################################################################################



per_batch_normalisation <- function(input.file = "./2_results/scNMT/RNAseq/sce_qc2.rds",
                                    output.file = "./2_results/scNMT/RNAseq/Batch_Correction/sce_qc2_PerBatchNormPCA.rds",
                                    output.fig.dir = './figures/scNMT/RNAseq/Batch_Correction/Normalisation/',
                                    per.batch = TRUE, # If false, it performs scuttle:logNormalise() for the combined batches but still performs runPCA on the batches separately after.
                                    force = FALSE){
  
  
  # Check existence and create dir if necessary
  
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  
  
  # Read in file
  cat("\nloading the input file:", input.file)
  sce_qc <- readRDS(input.file)
  unique(sce_qc$Batch)
  
  ## ------------- Multi-Batch Normalisation
  cat("\nperforming per-batch normalisation ...\n")
  cat("\nwith post normalisation PCA and umap calculations\n")
  sce.list <- list()
  
  if(per.batch == TRUE) {
    
    Batches <- unique(sce_qc$Batch)
    
    for (Batch in unique(sce_qc$Batch))
    {
      set.seed(405)
      sce_i <- sce_qc[,sce_qc$Batch == Batch]
      sce_i <- logNormCounts(sce_i)
      sce_i <- runPCA(sce_i, ncomponents = 4, name = 'PCA.batch')
      sce_i <- runUMAP(sce_i, ncomponents = 4, name = 'UMAP.batch')
      sce.list[[Batch]] <- sce_i
    }
  } else if (per.batch == FALSE){
  
    sce_qc <- logNormCounts(sce_qc)
    Batches <- unique(sce_qc$Batch)
    
    for (Batch in unique(sce_qc$Batch))
    {
      set.seed(405)
      sce_i <- sce_qc[,sce_qc$Batch == Batch]

      sce_i <- runPCA(sce_i, ncomponents = 4, name = 'PCA.batch')
      sce_i <- runUMAP(sce_i, ncomponents = 4, name = 'UMAP.batch')
      sce.list[[Batch]] <- sce_i
    }
  }  
  
  sce.batch <- Reduce('cbind', sce.list)
  
  # Plot norm factors used
  p1 <- as.data.frame(colData(sce.batch)) %>% ggplot(aes(x=total, y =sizeFactor, color = Batch, shape = Treatment)) + geom_point() + coord_trans(x="log2", y = "log2") + theme_light() + xlab("Library size")
  p2 <- as.data.frame(colData(sce.batch)) %>% ggplot(aes(x=total, y =sizeFactor, color = Treatment, shape = Batch)) + geom_point() + coord_trans(x="log2", y = "log2") + theme_light() + xlab("Library size")+
    scale_colour_manual(values = c("Unt" = "#00BA38", "AZA" = "#619CFF", "DAC" = "#F8766D"))
  plot_grid(p1, p2, labels = c('A', 'B'))
  if(per.batch==TRUE){
    ggsave(filename = paste0(output.fig.dir,"SizeFactor.MeanScale_vs_lib.size_perBatchTrue.png"), width = 18, height = 6)
  }else if(per.batch==FALSE){
    ggsave(filename = paste0(output.fig.dir, "SizeFactor.MeanScale_vs_lib.size_perBatchFALSE.png"), width = 18, height = 6)
  }
  
  # runPCA with batches together post normalisation
  sce.batch <- runPCA(sce.batch, ncomponents = 4, name = 'PCA.batch.united', exprs_values = "logcounts")
  
  # TODO variance explained by each PC
  cat("writing output file ...\n")
  message(output.file)
  saveRDS(sce.batch, output.file)
  
  cat("\nFinished normalisation script. Now run visualisaiton script.\n")
  
}



########################################################################################################
## per_batch_normalisation_CPM_NoScaling
########################################################################################################



per_batch_normalisation_CPM_NoScaling <- function(input.file = "./2_results/scNMT/RNAseq/sce_qc2.rds",
                                    output.file = "./2_results/scNMT/RNAseq/Batch_Correction/Normalisation/sce_qc2_PerBatchNormPCA_CPM_NoScaling.rds",
                                    output.fig.dir = './figures/scNMT/RNAseq/Batch_Correction/Normalisation/'){
  
  
  # Check existence and create dir if necessary
  
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  
  
  # Read in file
  cat("\nloading the input file:", input.file)
  sce_qc <- readRDS(input.file)
  unique(sce_qc$Batch)
  
  ## ------------- Multi-Batch Normalisation
  cat("\nperforming per-batch normalisation ...\n")
  cat("\nwith post normalisation PCA and umap calculations\n")
  sce.list <- list()
  
  # Perform CPM without scaling
  raw_counts <- assay(sce_qc)
  lib.size <- colSums(raw_counts)
  
  CPM_counts <- log2(t(t(raw_counts + 1)/(1e-06 * lib.size))) # Equivalent to log2(x+1/colsums(x) *1million) and t((t(raw_counts) +1)/lib.size)*1000000
  
  assay(sce_qc, "CPM_No_Scaling") <- CPM_counts
  
  Batches <- unique(sce_qc$Batch)
  
  for (Batch in unique(sce_qc$Batch))
  {
    set.seed(405)
    sce_i <- sce_qc[,sce_qc$Batch == Batch]
    
    sce_i <- runPCA(sce_i, ncomponents = 4, name = 'PCA.batch', exprs_values = "CPM_No_Scaling")
    sce_i <- runUMAP(sce_i, ncomponents = 4, name = 'UMAP.batch', exprs_values = "CPM_No_Scaling")
    sce.list[[Batch]] <- sce_i
  }
  
  
  sce.batch <- Reduce('cbind', sce.list)
  

  # runPCA with batches together post normalisation
  sce.batch <- runPCA(sce.batch, ncomponents = 4, name = 'PCA.batch.united', exprs_values = "CPM_No_Scaling")
  
  # TODO variance explained by each PC
  cat("writing output file ...\n")
  message(output.file)
  saveRDS(sce.batch, output.file)
  
  cat("\nFinished normalisation script. Now run visualisaiton script.\n")
  
}

########################################################################################################




per_batch_normalisation_Pooled_Deconvolution_DEVELOPMENT <- function(input.file = "./2_results/scNMT/RNAseq/sce_qc2.rds",
                                    output.file = "./2_results/scNMT/RNAseq/Batch_Correction_scran/sce_qc2_Norm-PoolDeconv.rds",
                                    min.size = 20,
                                    cluster.method = "igraph", # Default is igraph. Alternative is 'hclust' or 'NONE'. If none. min.size will be ignored (no need to change it).
                                    perBatch = FALSE, # But runPCA is still performed by batch
                                    force = FALSE){
  
  # This applies the pooled-based deconvolution approach from Lun et al 2016. Which was first written in the scran pacakge and migrated to the scuttle package.
  # As this relies on pooling groups of cells (default was pools of 100 cells) this is not written to allow splitting of batches for individual normalisation.
  # The PCA after normalisation are calculated on the batches independently to identify the most variable genes/TEs specific to their own batches.
  # A separate PCA is run on newly 'normalised' data to find the most variable genes - most likely this will pick up genes that explain differences in batches. Hence the requirement of additional batch correction methods...
  
  
  # Check existence and create dir if necessary
  
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  
  
  # Read in file
  cat("\nloading the input file:", input.file)
  sce_qc <- readRDS(input.file)
  unique(sce_qc$Batch)
  
  ## ------------- Multi-Batch Normalisation
  cat("\nperforming per-batch normalisation ...\n")
  cat("\nwith post normalisation PCA and umap calculations\n")
  sce.list <- list()
  
  
  # Computing the size factors
  # Using pooled with pre-clustering and scuttle/scran
  ## Perform Clustering of cells prior for pooled devonvolution approach.
  library(scran)
  preclusters <- quickCluster(sce_qc, min.size = 20, method = cluster.method) # Modified to 20 from default 50 to allow for 1/3 of a batch (roughly one treatment group) to be the minimal cluster size. Arguably this could slightly smaller...
  table(preclusters)
  # data.frame(Sample = colnames(sce_qc), Batch = colData(sce_qc)$Batch, Treatment = colData(sce_qc)$Treatment,  Cluster = preclusters) 
  # data.frame(Sample = colnames(sce_qc), Batch = colData(sce_qc)$Batch, Treatment = colData(sce_qc)$Treatment,  Cluster = preclusters) %>% filter(Cluster ==2)
  as.data.frame(colData(sce_qc_scuttle)) %>% ggplot(aes(x=total, y =sizeFactor, color = Batch, shape = Treatment)) + geom_point()
  ## Compute N
  sce_qc_scuttle <- computePooledFactors(sce_qc) # No Clustering because sample size so small...?
  sce_qc_scuttle$sizeFactor.deconv2 <- scuttle::pooledSizeFactors(sce_qc_scuttle, scaling=sce_qc_scuttle$sizeFactor)
  
  sce_qc_scuttle2 <- computePooledFactors(sce_qc, sizes = 15, clusters=preclusters) # this equivalent to 'computeSumFactors(sce_qc, sizes = 15, clusters = preclusters)'
  head(sizeFactors(sce_qc2_scuttle))
  sce_qc_scuttle2$Clusters <- preclusters
  as.data.frame(colData(sce_qc_scuttle2)) %>% ggplot(aes(x=total, y =sizeFactor, color = Batch, shape = Treatment)) + geom_point()
  as.data.frame(colData(sce_qc_scuttle2)) %>% ggplot(aes(x=total, y =sizeFactor, color = Clusters, shape = Treatment)) + geom_point()
  
  ## Use Calculated Norm Factors to log2 normalise
  sce_qc_scran <- scater::logNormCounts(sce_qc_scran)

  

  
  Batches <- unique(sce_qc$Batch)
  for (Batch in unique(sce_qc$Batch))
  {
    set.seed(405)
    sce_i <- sce_qc[,sce_qc$Batch == Batch]
    sce_i_clusters <- quickCluster(sce_i)
    sce_i <- computeSumFactors(sce_i, clusters=sce_i_clusters)
    
    lib.sf <- librarySizeFactors(sce_i)
    summary(lib.sf)
    edgeR::calcNormFactors(sce_i, method="RLE")
    
    dgeSeq = edgeR::DGEList(counts=assay(sce_i, 'counts'))
    dgeSeq = edgeR::calcNormFactors(dgeSeq, method="RLE")
    sce_i <- computeSumFactors(sce_i)
    sizeFactors(sce_i)
    
    plot(x= dgeSeq$samples$lib.size, y = dgeSeq$samples$norm.factors)
    plot(x= dgeSeq$samples$lib.size, y = sizeFactors(sce_i))
    
    assay(sce_i, 'logcounts2', withDimnames=F) = edgeR::cpm(dgeSeq, log=T)
    dec2 <- modelGeneVar(sce_i, assay.type ='logcounts2')
    plot(dec2$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
    curve(metadata(dec2)$trend(x), col="blue", add=TRUE)
    
    sce_i <- logNormCounts(sce_i)
    
    dec <- modelGeneVar(sce_i, assay.type ='logcounts')
    plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
    curve(metadata(dec)$trend(x), col="blue", add=TRUE) 
    
    par(mfrow=c(1,2))
    plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance")
    curve(metadata(dec)$trend(x), col="blue", add=TRUE) 
    plot(dec2$mean, dec2$total, xlab="Mean log-expression", ylab="Variance")
    curve(metadata(dec2)$trend(x), col="blue", add=TRUE)
    
    sce_i <- runPCA(sce_i, ncomponents = 4, name = 'PCA.batch')
    sce_i <- runUMAP(sce_i, ncomponents = 4, name = 'UMAP.batch')
    sce.list[[Batch]] <- sce_i
  }
  
  sce.batch <- Reduce('cbind', sce.list)
  
  # TODO variance explained by each PC
  cat("writing output file ...\n")
  message(output.file)
  saveRDS(sce.batch, output.file)
  
  cat("\nFinished normalisation script. Now run visualisaiton script.\n")
  
}