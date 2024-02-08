##########################################################################
# Purpose: R Script to perform Seurat4 integration/batch correction
# Output: Seurat corrected RDS and figures
#
# Date: 25.August.22
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) Requires pre-normalised (but QCd) data in a SingleCellExperiment Object.
## 2) Requires 'Batch' column in metadata of SingleCellExperiment object.
## 3) This script automatically creates figures - calls in secondary helper scripts.
# To do: 
## 1) Fix SCE-> seurat-> SCE conversion: Seurat replaces underscores ('_'), with dashes ('-'). Automate this 're-correction' in updated sce. Can get initial list with: row.names(sce)[grepl("_",row.names(sce))]
##
##########################################################################
# Relevant online material:
## https://satijalab.org/seurat/articles/integration_introduction.html
## https://satijalab.org/seurat/articles/sctransform_v2_vignette.html
## https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8 (original method described)
# http://bioconductor.org/books/3.13/OSCA.advanced/integrating-with-protein-abundance.html#integrating-with-protein-abundance
## Consider info on imputation for seurat method: https://satijalab.org/seurat/archive/v3.1/atacseq_integration_vignette
## https://pubmed.ncbi.nlm.nih.gov/31178118/
## https://biocellgen-public.svi.edu.au/mig_2019_scrnaseq-workshop/comparing-and-combining-scrna-seq-datasets.html#anchor-integration-seurat
## Preserving zeros: https://rdrr.io/github/satijalab/seurat-wrappers/man/RunALRA.html
## Comparison of multiple imputation methods: https://rnabioco.github.io/cellar/previous/2019/docs/8_imputation.html
## http://bioconductor.org/books/release/OSCA/
## https://www.singlecellcourse.org/scrna-seq-dataset-integration.html#seurat-v3-3-vs-5-10k-pbmc
## chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://trace.tennessee.edu/cgi/viewcontent.cgi?article=8909&context=utk_graddiss
## https://github.com/satijalab/seurat/issues/1732
## Conversion on non-zero controversy: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02601-5#Sec6
## https://swaruplab.bio.uci.edu/tutorial/integration/integration_tutorial.html#sctransform
## https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html
## http://bioconductor.org/books/3.12/OSCA/cell-cycle-assignment.html
## https://hbctraining.github.io/scRNA-seq/lessons/06_SC_SCT_and_integration.html
## https://htmlpreview.github.io/?https://github.com/satijalab/sctransform/blob/supp_html/supplement/variance_stabilizing_transformation.html
##########################################################################
# Load packages
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(ggplot2, 
               knitr, 
               SingleCellExperiment,
               mixOmics, 
               ggplot2,
               batchelor,
               rlang,
               scater,
               Seurat,
               gsubfn)

source("./scripts/converted_from_HPC/viz_batch_normalisation.R") # For custom plotting

batch_correction_seurat4 <- function(input.file = "2_results/scNMT/RNAseq/sce_qc2.rds", 
                                     output.file = "2_results/scNMT/RNAseq/Batch_Correction/Seurat/Seurat4.rds",
                                     output.fig.dir = 'figures/scNMT/RNAseq/Batch_Correction/',
                                     assay.name = "Seurat", # This is the assay name in the output RDS.
                                     batch.column = 'Batch',
                                     batch.order = NA, # This forces plots to display in the specified order for the specified batch.column c("b320","b620", "b322")
                                     seurat.norm = "SCTransform", # Can be either 'SCTransform' or 'vst'.
                                     seurat.SCTransform.vst.flavor = "v2", # NOT CODED YET. Only used if 'seurat.norm = SCTransform'. 
                                     nhvgs = 5000, # Number of highly variable genes used for anchoring detection and integration (retained in output batch corrected rds file)
                                     k.weight = 50, # seurat default is 100. Adjusted for lower throughput scNMT data (96 well plates + QC requirements).
                                     set.the.seed = TRUE,
                                     force = FALSE){
 
  ## Set seed (advisable for PCA plotting etc. reproducibility)
  if(set.the.seed == TRUE){
    set.seed(123)
  }
  
  ## Create output dir (if not already present)
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(output.fig.dir), recursive = TRUE, showWarnings = FALSE)
  
  # Read in sce
  sce <- readRDS(input.file)
  # Convert to Seurat object
  
  sce.seurat <- CreateSeuratObject(counts = assay(sce,"counts"), 
                                   project = "seurat.temp", min.cells = 0,
                                   meta.data = as.data.frame(colData(sce)))
  ## Get list of any features that have underscores replaced (for later 're-conversion' back to SCE object)
  SCE_to_Seurat_underscores <- row.names(sce)[grepl("_",row.names(sce))]
  SCE_to_Seurat_underscores_removed <- gsub("_","-",SCE_to_Seurat_underscores) # use gsub to match seurat object for ID to convert later. i.e. if(x %in% SCE_to_Seurat_underscores_removed){convert to SCE_to_Seurat_underscores}
  
  ## by SCTransform normalisation
  ########################################################################################################
  
  if(seurat.norm == "SCTransform"){
  
    # split the dataset into a list of separate seurat objects by 'Batch'
    ifnb.list <- SplitObject(sce.seurat, split.by = "Batch")
    
    # Identify variable features for each dataset independently
    ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform)
    
    # select features that are repeatedly variable across datasets for integration
    features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = nhvgs)
    ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
  
    # this command creates an 'integrated' data assay
    immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                            anchor.features = features)
    immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT",  k.weight = k.weight)

    # Run the standard workflow for visualization and clustering
    immune.combined.sct <- RunPCA(immune.combined.sct, npcs = 30, verbose = FALSE, seed.use = 123, features =features)
    immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30, seed.use = 123)
    
    # Visualization
    ## Force order of groups plotted
    if(is.na(batch.order) != TRUE){ 
      immune.combined.sct[[batch.column]] <- factor(unlist(as.data.frame(t(immune.combined.sct[[batch.column]]))), levels = batch.order) # Slightly more complex approach required to modify factor levels using a variable colname with a seurat object.
    }
    ## Combined plot (Seurat DimPlot)
    ### UMAPs
    p1 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "Batch", pt.size = 2.5)
    p2 <- DimPlot(immune.combined.sct, reduction = "umap", group.by = "Treatment", pt.size = 2.5)
    p3 <- DimPlot(immune.combined.sct, reduction = "umap", split.by = "Batch", group.by = "Treatment", pt.size = 2.5)
    ### PCAs
    p4 <- DimPlot(immune.combined.sct, reduction = "pca", group.by = "Batch", pt.size = 2.5)
    p5 <- DimPlot(immune.combined.sct, reduction = "pca", group.by = "Treatment", pt.size = 2.5)
    p6 <- DimPlot(immune.combined.sct, reduction = "pca", split.by = "Batch", group.by = "Treatment", pt.size = 2.5)
    ### Combine plots
    p_row1 <- cowplot::plot_grid(p1, p2, ncol = 2)
    p_row3 <- cowplot::plot_grid(p4, p5, ncol = 2)
    cowplot::plot_grid(p_row1, 
                     p3,
                     p_row3,
                     p6,
                     ncol = 1, nrow = 4)
    
    ### Save plot
    output.fig.name = paste0("sce_qc2_Seurat4_",seurat.norm,"_",nhvgs,"_batch_vs_tx.png")
    filename <- paste0(output.fig.dir, "/", output.fig.name)
    cat("\nwriting: ", filename, "\n")
    ggsave(
      filename =  filename,
      width = 12,
      height = 20)
  
    
    ## Customised ggplot with side distribution plots.
    # source("./scripts/converted_from_HPC/viz_batch_normalisation.R")
    ### Convert to SCE object
    pbmc.sce <- as.SingleCellExperiment(immune.combined.sct)
    #pbmc.sce2 <- runPCA(pbmc.sce, exprs_values = "logcounts", name = "PCA_scater") # Default uses 'top 500 hvgs' identified by scater package.
    
    ### PCA (Seurat values)
    shape = "Batch"
    size = 3
    width = 12
    height = 8
    
    plot_reducedDim_Norm(
      sce = pbmc.sce,
      reducedDim = "PCA",
      comps = c(1, 2),
      colBy = 'Treatment',
      cols = c("Unt" = "#00BA38", "AZA" =  "#619CFF", "DAC"  = "#F8766D"),
      axis.title = 'PC',
      size = size,
      shape = "Batch"
    )
    
    p_to_save <- grid.arrange(Plots_saved$xdensity, Plots_saved$blankPlot, Plots_saved$scatter, Plots_saved$ydensity, 
                              ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
    
    #output.fig = 'figures/scNMT/RNAseq/Batch_Correction/' # output folder for figures
    output.fig.name = paste0("sce_qc2_Seurat4_",seurat.norm,"_",nhvgs,"_PCA.png")
    
    filename <- here(paste0(output.fig.dir, "/", output.fig.name))
    cat("\nwriting: ", filename, "\n")
    ggsave(plot  = p_to_save, 
           filename =  filename,
           width = width,
           height = height)
    
    
    ## UMAP (Seurat values)
    
    plot_reducedDim_Norm_UMAP(
      sce = pbmc.sce,
      reducedDim = "UMAP",
      comps = c(1, 2),
      colBy = 'Treatment',
      cols = c("Unt" = "green", "AZA" =  "blue", "DAC"  = "red"),
      axis.title = 'UMAP',
      size = size,
      shape = shape
    )
    
    p_to_save <- grid.arrange(Plots_saved$xdensity, Plots_saved$blankPlot, Plots_saved$scatter, Plots_saved$ydensity, 
                              ncol=2, nrow=2, widths=c(4, 1.4), heights=c(1.4, 4))
    
    #output.fig = 'figures/scNMT/RNAseq/Batch_Correction/' # output folder for figures
    output.fig.name = paste0("sce_qc2_Seurat4_",seurat.norm,"_",nhvgs,"_UMAP.png")
    
    #filenam <- here(sprintf('%s/colony-PerBatch_QC2_PCA.png', output.fig))
    filename <- here(paste0(output.fig.dir, "/", output.fig.name))
    
    cat("\nwriting: ", filename, "\n")
    ggsave(plot  = p_to_save, 
           filename =  filename,
           width = width,
           height = height)
    
  }
  # Save rds file (for now with Seurat object)
  cat("Saving seurat object in rds file: ",output.file)
  seurat.corrected.sct <- immune.combined.sct # Changing object name. Modify this above to save adjusting it here.
  saveRDS(seurat.corrected.sct, file = output.file)
}    


seurat_update_with_NAs <- function(input.seurat.rds = "2_results/scNMT/RNAseq/Batch_Correction/Seurat/Seurat4.rds",
                                   output.file = "2_results/scNMT/RNAseq/Batch_Correction/Seurat/Seurat4_NAs.rds",
                                   assay.with.missing.values = "RNA",
                                   assay.to.add.missing.values = "integrated",
                                   new.assay.name = "integrated_NAs"
                                   ){

  ## Create output dir (if not already present)
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)

  # Read in file
  seurat.orignal <- readRDS(input.seurat.rds)

  # Extract assay to add missing values to (integrated)
  assay.needs.NAs <- GetAssayData(seurat.orignal[[assay.to.add.missing.values]]) %>% as.matrix()

  # Extract assay with missing values and get index of genes with missing values/zeros
  ## Assay with missing values
  assay.with.missing.values <- as.matrix(GetAssayData(seurat.orignal[[assay.with.missing.values]])) %>% as.matrix()
  ## Filter for only genes present in the second assay (with one that needs NAs added)
  assay.with.missing.values_subset <- assay.with.missing.values[rownames(assay.needs.NAs),]

  ## Make sure row and columns are in the same order!
  cat("\nChecking if rownames (features) are in the same order:\n")
  identical(rownames(assay.with.missing.values_subset),rownames(assay.needs.NAs))
  cat("\nChecking if colnames (samples) are in the same order:\n")
  identical(colnames(assay.with.missing.values_subset),colnames(assay.needs.NAs))
  cat("\ndelete output if either were false!\n") # Add stop if above.

  ## Get missing value index
  missing_value_index <- which(assay.with.missing.values_subset == 0, arr.ind = T)

  cat("\nidentified a total of ", nrow(missing_value_index), "missing values to add.\n")
  Percent.missing <- (nrow(missing_value_index)/ (nrow(assay.needs.NAs) *ncol(assay.needs.NAs)) *100)
  cat("\nThis accounts for", Percent.missing,
      "% in the updated assay.\n")
  ## Compare to NAs present in QCd data pre-normalisation and integration
  missing_value_index_2 <- which(assay.with.missing.values == 0, arr.ind = T)
  cat("\nFor comparison, identified a total of ", nrow(missing_value_index2), "missing in raw data.\n")
  Percent.missing_2 <- (nrow(missing_value_index2)/ (nrow(assay.with.missing.values) *ncol(assay.with.missing.values)) *100)
  cat("\nThis accounts for", Percent.missing_2,
      "% in the pre-norm/corrected data.\n")

  # Update matrix with NAs
  assay.updated <- assay.needs.NAs
  assay.updated[missing_value_index] <- NA
  ## Check correct number of NAs we created
  nrow(missing_value_index) == sum(is.na(assay.updated))

  # Add updated assay to seurat object
  updatedAssay.seurat <- CreateAssayObject(counts = assay.updated)
  seurat.orignal[[new.assay.name]] <- updatedAssay.seurat

  # Save updated seurat object
  cat("Saving updated seurat object in rds file: ",output.file)
  saveRDS(seurat.orignal, file = output.file)
}


