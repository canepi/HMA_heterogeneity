##########################################################################
# Purpose: R Script to perform scMerge correction
# Output: scMerge corrected sce object saved as an RDS file
#
# Date: 06.Mar.23
# Version: v.0.2.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: This was originally developed as part of a package build running on unix/HPC by Al, which would allow shell style submission.
## 1) added assay.type for variable input/selection - defaults still to 'logcounts'.
# To do: Enable both parallelisation with svd - should be quicker but yields slightly different (albeit very similar) results. The 'normal' version can also be run with parallelisation too. Which is already enabled here.
##########################################################################
if(!require(pacman)){install.packages("pacman");require(pacman)}
# pacman::p_unload(all) # This is important as I found clashes between packages when running 'runPCA'.

pacman::p_load(ggplot2, 
               knitr, 
               SingleCellExperiment,
               mixOmics, 
               ggplot2,
               scMerge,
               rlang,
               scater)


batch_correction_scMerge <- function(
    input.file = '2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect.rds',
    output.file = "2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_3C.rds",
    assay.name = "scMerge3C", # Name will concatenated with '_celltype' for correction using cell type info. Recommend add number of clusters provided 'XC',
    assay.input = "logcounts", # assay name input
    Cluster.number = 3, # Number of 'kmeans' clusters. Suggestion would be number of treatments or cell types.
    Batch.number = 3, # Number of batches being measures. The cluster.number and batch.number will be used in combination to form the 'kmeansk' variable for scMerge().
    force = FALSE,
    parallel = TRUE,
    fast.svd = "Y"){ # If you select 'YES' then it will use svd option. Which should run quicker, but yield slightly different (albeit similar) results to 'normal' approach. Only use with extremely large datasets (>3,000 cells).


if (force | !file.exists(output.file)){
#  input.file = '2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect.rds'
#  output.file = "2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_3C.rds"
#  assay.name = "scMerge3C"
#  Cluster.number = 3
#  Batch.number = 3
#  force = FALSE
#  parallel = TRUE
#  fast.svd = "Y"

# Create dir
dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)

# Generate kmeansK from cluster.number and batch.number provided.
KmeansK_var <- rep(Cluster.number, Batch.number)



cat("\nloading the input data:\n", input.file, "\n")
sce <- readRDS(input.file)
# https://sydneybiox.github.io/scMerge/articles/scMerge.html#selectnclibrary(scMerge)
## SEG list in ensemblGene ID
data("segList_ensemblGeneID", package = "scMerge")
## SEG list for human scRNA-Seq data
head(segList_ensemblGeneID$human$human_scSEG)



######################################################################################################################################################################
if(parallel == TRUE && fast.svd %in% c("NO","N")){
  message("Running in parallel with 'normal' scMerge algorithm (not svd)")

# Register parallelisation options to work across platforms (windows and unix). requires SnowParam.
# For more information on parallelisaiton and svd method to speed up, see: https://sydneybiox.github.io/scMerge/articles/scMerge.html#achieving-fast-and-memory-efficient-computation
  library(BiocParallel)
  cores_to_use <- parallel::detectCores()-1
  message("\nRunning scMerge in Parallel.\n", cores_to_use," cores will be used, which is 1 less core than the total available. " )

## without inclusion of cell-type info & 'original'/slower approach (non-svd)
  cat("\nFirst, Performing scMerge WITHOUT celltype info and saving assay to SCE...\n")
  assay.name <- assay.name
  PCA.name <- paste0('PCA_', assay.name)
  
  sce_Merge <- scMerge(
    sce_combine = sce,
    ctl = (segList_ensemblGeneID$human$human_scSEG),
    exprs = assay.input, # default is 'logcounts'
    batch_name = "Batch",
    # cell_type = sce$Treatment,
    # cell_type_match = TRUE,
    kmeansK = KmeansK_var,
    assay_name = assay.name,
    BPPARAM = SnowParam(workers = cores_to_use)
    )

  
## with inclusion of cell-type info
  cat("\nPerforming scMerge USING celltype info and saving assay to SCE...\n")
  
  cta <- paste0(assay.name, "_celltype")
  PCA.name.cta <- paste0(PCA.name, "_celltype")
  
  sce_Merge <- scMerge(
    sce_combine = sce_Merge,
    ctl = (segList_ensemblGeneID$human$human_scSEG),
    exprs = "logcounts",
    batch_name = "Batch",
    cell_type = sce$Treatment,
    cell_type_match = TRUE,
    kmeansK = KmeansK_var,
    assay_name = cta,
    BPPARAM = SnowParam(workers = cores_to_use)
    )
  
  
  # Running PCA on scMerged results

  cat("\nperforming PCA on scMerged data ...\n")
  sce_Merge <- runPCA(sce_Merge, exprs_values = assay.name, name = PCA.name)
  sce_Merge <- runPCA(sce_Merge, exprs_values = cta,        name = PCA.name.cta)
  
  cat("\nsaving the output SCE...\n",output.file)
  
  saveRDS(sce_Merge, file = output.file)
  }

##########################################################################################################
if(parallel == TRUE && fast.svd %in% c("YES","Y")){
  message("Running in parallel with 'faster svd approach'.")

  library(BiocSingular)
  
## without inclusion of cell-type info & 'faster' svd approach
  cat("\nFirst, Performing scMerge WITHOUT celltype info and saving assay to SCE...\n")
  
  assay.name <- assay.name
  PCA.name <- paste0('PCA_', assay.name)

  sce_Merge <- scMerge(
    sce_combine = sce,
    ctl = (segList_ensemblGeneID$human$human_scSEG),
    exprs = assay.input, # default is 'logcounts'
    batch_name = "Batch",
    # cell_type = sce$Treatment,
    # cell_type_match = TRUE,
    kmeansK = KmeansK_var,
    assay_name = assay.name,
    BSPARAM = ExactParam()
    )


## with inclusion of cell-type info
  cat("\nPerforming scMerge USING celltype info and saving assay to SCE...\n")

  cta <- paste0(assay.name, "_celltype")
  PCA.name.cta <- paste0(PCA.name, "_celltype")
  
  sce_Merge <- scMerge(
    sce_combine = sce_Merge,
    ctl = (segList_ensemblGeneID$human$human_scSEG),
    exprs = "logcounts",
    batch_name = "Batch",
    cell_type = sce$Treatment,
    cell_type_match = TRUE,
    kmeansK = KmeansK_var,
    assay_name = cta,
    BSPARAM = ExactParam()
    )
  
  
# Running PCA on scMerged results
    # Need to remove clashing package with 'runPCA' (found one of the packages that enable the parallelisation clashes with scater runPCA/BiocSingular - so probably is BiocParallel)
    # pacman::p_unload(BiocParallel)
  cat("\nperforming PCA on scMerged data ...\n")
  sce_Merge <- runPCA(sce_Merge, exprs_values = assay.name, name = PCA.name)
  sce_Merge <- runPCA(sce_Merge, exprs_values = cta,        name = PCA.name.cta)
  
  cat("\nsaving the output SCE...\n")
  
  saveRDS(sce_Merge, file = output.file)
  }
##################################################################################################################################################################################################
# Running without parallelism
if(parallel == FALSE){
  message("\nRunning scMerge on a single core. This might take some time with large datasets.,")

  cat("\nFirst, Performing scMerge WITHOUT celltype info and saving assay to SCE...\n")
  assay.name <- assay.name
  PCA.name <- paste0('PCA_', assay.name)

  sce_Merge <- scMerge(
    sce_combine = sce,
    ctl = (segList_ensemblGeneID$human$human_scSEG),
    exprs = assay.input, # default is 'logcounts'
    batch_name = "Batch",
    # cell_type = sce$Treatment,
    # cell_type_match = TRUE,
    kmeansK = KmeansK_var,
    assay_name = assay.name
    )
  
  ## with inclusion of cell-type info
  cat("\nPerforming scMerge USING celltype info and saving assay to SCE...\n")
  
  cta <- paste0(assay.name, "_celltype")
  PCA.name.cta <- paste0(PCA.name, "_celltype")
  
  sce_Merge <- scMerge(
    sce_combine = sce_Merge,
    ctl = (segList_ensemblGeneID$human$human_scSEG),
    exprs = assay.input, # default is 'logcounts'
    batch_name = "Batch",
    cell_type = sce$Treatment,
    cell_type_match = TRUE,
    kmeansK = KmeansK_var,
    assay_name = cta
    )
  
  
  # Running PCA on scMerged results
  cat("\nperforming PCA on scMerged data ...\n")
  sce_Merge <- runPCA(sce_Merge, exprs_values = assay.name, name = PCA.name)
  sce_Merge <- runPCA(sce_Merge, exprs_values = cta,        name = PCA.name.cta)
  
  cat("\nsaving the output SCE...\n")
  
  saveRDS(sce_Merge, file = output.file)
  }
} else {message("Output file existance = ", file.exists(output.file),"\nChange 'force =TRUE' to overwrite.")}
}

