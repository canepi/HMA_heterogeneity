##########################################################################
# Purpose: R Script to perform NIPALS on acc/met assays
# Output: RDS files containing the PCA calculations and associated PCA plots (png)
#
# Date: 29.03.23
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) 
# To do: 
## 1) Option to select contexts/assays to calculate
## 2) Allow automatic imputation when detecting empty PCA calcuated. And record which ones were imputed!
# Websites:
##
##########################################################################
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(ggplot2, dplyr, knitr, SingleCellExperiment, rlang, Seurat, nipals, MultiAssayExperiment, doParallel)


calc_pca_met_acc <- function(input.file = "./2_results/scNMT/NMT/NMT_MAE.rds",
                             output.folder = "./2_results/scNMT/NMT/pca",
                             n_hvrs = 5000,
                             batch = NA,
                             pctcells = 50,
                             ncomps = 3,
                             fitted = TRUE,
                             run.parallel = FALSE,
                             assays = NA, # Defaults to all. Specifying assays can help reduce time taken due to imputation (if required).
                             allow.imputation = FALSE, # Warning. This can take a while for a large datasets and if many assays need to be imputed...
                             set.the.seed = TRUE){


  # Set seed (for reproducibility of PCA plots)
  if(set.the.seed == TRUE){set.seed(123)}

  cat("\nusing data path: ", input.file)
  cat("\nusing output path: ", output.folder)
  cat("\nusing nhvr: ", n_hvrs)
  pct_cells_detected <- pctcells
  cat("\nusing pct_cells_detected: ", pct_cells_detected)
  cat("\nusing ncomps: ", ncomps)
  cat("\nusing fitted: ", fitted)


  cat("\n\nimporting the MAE ...")
  HMA.MAE <- readRDS(input.file)

  
  # Filtering for selected batch (if required)
  if (any(is.na(batch)) == FALSE) { # subset if batch is provided. If NA. Keeps all samples
    HMA.MAE <- HMA.MAE[,HMA.MAE$Batch %in% batch]
  }
  cat("\nKept batches:", unique(HMA.MAE$Batch))
  
  assay.names <- names(assays(HMA.MAE))

  ## perform PCA on all assays
  cat("\nCalculating PCAs")
  
  if(run.parallel == TRUE){
    cat("\nRun parallel =", run.parallel)
    cat("\ninitialising doParallel ...\n")
  
    N_CORES <- parallel::detectCores() # Detect how many cores are available
    my.cl <- parallel::makeCluster((N_CORES-2), type = "PSOCK") # Make cluster using 'socket' which works across unix and windows! Using 1 core less than totally available to allow for other work to be done while waiting.
    doParallel::registerDoParallel(cl = my.cl) # Register parallel options
  }
  
  out.dir <- file.path(output.folder, sprintf("Batch_%s_pctCells_%s_nhvrs_%s_ncomps_%s", paste0(batch, collapse = '_'), pct_cells_detected, n_hvrs, ncomps))
  fnames.out.pca <- file.path(out.dir, paste0(assay.names,'.rds'))
  fnames.out.fitted <- file.path(out.dir, 'fitted', paste0(assay.names,'.rds'))
  
  #For practice # i = 1:length(assay.names)
  #i = 1
  cat("\nCalculating NIPALS for:", as.vector(rbind('\n', assay.names)))
  
  foreach(i = 1:length(assay.names), .inorder=FALSE) %do% {
  
    assay.name <- assay.names[i]
    cat("\nPerforming NIPALS on: ... ", assay.name)
  
    fname.out.pca <- fnames.out.pca[i]
    fname.out.fitted <- fnames.out.fitted[i]
  
    if (!dir.exists(dirname(fname.out.pca)))
      dir.create(dirname(fname.out.pca), recursive = TRUE)
    if (!dir.exists(dirname(fname.out.fitted)))
      dir.create(dirname(fname.out.fitted), recursive = TRUE)
    
    assay.value <- assay(HMA.MAE, assay.name)
    
    # keep features detecting at least 50% of cells
    keep_features <- rowSums(!is.na(assay.value))/ncol(assay.value) >= pct_cells_detected/100
    assay.value <- assay.value[keep_features,]
    # keep top 5000 Highly Variable Regions
    hvrs <- rank(-rowVars(assay.value, na.rm = TRUE)) <= n_hvrs
    assay.value <- assay.value[hvrs,]
    nipals.res <- nipals::nipals(t(assay.value), ncomp = ncomps, fitted = fitted)
    
    if(fitted){
      saveRDS(nipals.res$fitted, file = fname.out.fitted)
      nipals.res$fitted <- NULL
      }
    
    saveRDS(nipals.res, file = fname.out.pca)
    
  }
  cat("\nFinished all PCA calculations with nipals.")
}
