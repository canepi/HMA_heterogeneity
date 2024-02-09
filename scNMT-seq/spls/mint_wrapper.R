##########################################################################
# Purpose: R Script to create MAE from CpG and GpC pre-computed rates
# Output: RDS file containing the results from the MixOmics analysis
#
# Date: 22.June.23
# Version: v.0.5.0
# Written by: Sean Burnard, adapted from Al J Abadi's original code.
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) Adapted from Al's original HPC script
## 2) Added options to select 'Experiments/blocks', assay (for all blocks), rna experiment and assay.
## 3) Added out.dir option.
## 4) Added keepY as vector (variable number of features per comp can be set).
## 5) Added keepXs as vector (variable number of features per comp can be set). But same applied to all X block.
# To do: 
## 1) 
# Websites:
##
##########################################################################
# Load packages
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(ggplot2, 
               dplyr,
               knitr, 
               SingleCellExperiment,
               data.table,
               #HDF5Array,
               SummarizedExperiment,
               MultiAssayExperiment,
               #ggplot2,
               #batchelor,
               rlang,
               scater,
               gsubfn, 
               mixOmics)


####################################################

mint_wrapper <- function(fun, # Should work with most mixOmics functions - (mint.)block.(s)pls(da) etc.
                         mae.path, 
                         out.dir = NA, 
                         assay.blocks.to.measure = NA, # Default will use all available blocks/assay. 
                         se.blocks.assay = NA, # Filters specific assay within all SingleExperiments used as blocks. Defaults to the first assay if multiple are available and none specified.
                         rna.experiment = 'rna',
                         rna.assay, # Specify which assay within the RNA experiment.
                         mode = 'canonical', 
                         pct_cells_detected,
                         n_hvrs, 
                         batch = NA, 
                         ncomp, 
                         keepXs = 50, 
                         force = FALSE, 
                         keepY=50) { # This can be a vector with different numbers per comp. If keepY != legnth(ncomp) the first value of keepY will be set for all comps.
  
  # perform mint.block.spls(da) and save the results
  DA <- grepl('da$', fun)
  if (DA)
  {
    mode <- NA
  }
  
  # Setting up output directory
  if (any(is.na(batch)) == FALSE) { 
    batches <- paste0(batch, collapse = '_')
  } else {batches = "NA"}
  ## If no out.dir given, use MAE dir
  if (is.na(out.dir) == TRUE) { 
    out.dir <- dirname(mae.path)
  }
  fname.out <- sprintf( paste(out.dir,"%s_%s_rna_%s%s_pctCells_%s_nhvrs_%s_ncomps_%s%s.rds", sep= ''), 
                        fun, batches, rna.assay, ifelse(is.na(se.blocks.assay), '', paste0('_SEblocksassay_', se.blocks.assay)), 
                        pct_cells_detected, n_hvrs, ncomp, sprintf('_mode_%s',mode))
  
  
  cat("\nusing output file: ", fname.out, '\n')
  
  if (force | !file.exists(fname.out))
  {
    #source("src/helpers/library.R")
    dir.create(dirname(fname.out), recursive = TRUE, showWarnings = FALSE)
    
    L1 <- c("\n\nimporting the MAE ...\n", mae.path)
    cat(L1)
    HMA.MAE <- readRDS(mae.path)
    
    # If batches specified, subset
    if (any(is.na(batch)) != TRUE) { # subset if batch is provided. If NA. Keeps all samples
      stopifnot(all(batch %in%  unique(HMA.MAE$Batch)))
      HMA.MAE <-   HMA.MAE[,HMA.MAE$Batch %in% batch]
    }
    cat("Kept batches:", as.character(unique(HMA.MAE$Batch)))
    
    # Extracting RNA assay
    rna <- HMA.MAE[[rna.experiment]] 
    rna <- t(assay(rna, rna.assay))
    # If experiments/blocks are specified
    if (any(is.na(assay.blocks.to.measure)) != TRUE) { # subset if batch is provided. If NA. Keeps all samples
      stopifnot(all(assay.blocks.to.measure %in%  names(assays(HMA.MAE))))
      HMA.MAE <- subsetByAssay(HMA.MAE, assay.blocks.to.measure)
    }
    assay.names <- names(assays(HMA.MAE))
    assay.names <- assay.names[!grepl(rna.experiment, assay.names)]
    names(assay.names) <- assay.names # for lapply
    
    L2 <- c("\nApplying MixOmics",fun,"using the following experiments as 'blocks':\n", as.vector(rbind('\n', assay.names)))
    cat(L2)
    L3 <- c("\nwith the RNA assay:", rna.assay)
    cat(L3)
    
    X <- lapply(assay.names, function(assay.name) {
      exp <- HMA.MAE[[assay.name]]
      if (is(exp, "SummarizedExperiment"))
      {
        exp <- as.matrix(assay(exp, se.blocks.assay))
      }
      mat <- as.matrix(t(exp))
    })
    
    X <- lapply(X, function(x) {
      keep <- colSums(!is.na(x))/nrow(x)*100 >= pct_cells_detected
      x[,keep]
    })
    
    X <- lapply(X, function(x) {
      mat <- x
      if (n_hvrs > 0 && dim(mat)[2] > n_hvrs)
      {
        hvrs <- rank(-colVars(mat, na.rm = TRUE)) <= n_hvrs
        mat <- mat[,hvrs]
      } else
      {
        mat
      }
    })
    
    X <- X[sapply(X, is.matrix)]
    ## keep blocks with at least 200 features
    drop_blocks <- names(X)[sapply(X, function(mat) dim(mat)[2] < 200)]
    if (length(drop_blocks))
      cat("\ndropping blocks with less than 200 features after filtering:", paste0(drop_blocks, collapse = ", "))
    X <- X[!names(X) %in% drop_blocks]
    
    DA <- grepl('da$', fun)
    
    if (DA)
    {
      Y <- HMA.MAE$Treatment
      rna <- rna[rownames(X[[1]]),] # Make sure rownames (samples) order of Expr matrix matches order in NOMeSeq data (assuming all other blocks are the same)
      X$rna <- as.matrix(rna)
    } else
    {
      rna <- rna[rownames(X[[1]]),] # Make sure rownames (samples) order of Expr matrix matches order in NOMeSeq data (assuming all other blocks are the same)
      Y <- as.matrix(rna)
      
      if(length(keepY) != ncomp) { # Allows flexibility of setting keepY. If ncomp doesn't match the length of keepY provided, then the first value of keepY will be set for all comps.
        cat("Automatically set keepY to:", rep(keepY, ncomp))
        keepY <- rep(keepY[1], ncomp)
      }
    }
    
    design <- matrix(0, nrow = length(X), ncol = length(X), dimnames = list(names(X), names(X)))
    
    # Setting KeepX 
    if(length(keepXs) != ncomp) { # Allows flexibility of setting keepY. If ncomp doesn't match the length of keepY provided, then the first value of keepY will be set for all comps.
      cat("Automatically setting keepX to:", rep(min(keepXs), ncomp))
      keepX <- lapply(X,  {
        keepXs <- min(keepXs, dim(x)[2])
        rep(keepXs, ncomp)
      })
    } else if(length(keepXs) == ncomp){
      keepX <- lapply(X, function(x){ # Assigns specified features per (X) comp, if length is equivalent to selected ncomps.
        keepXs})
    }

  
    ## get function symbol for call
    func <- get(fun)
    if (grepl('mint', fun))
    {
      # DA analysis uses regression mode only
      if (DA){
        L4 <- c("\nPerforming mint.block.splsda")
        cat(L4)
        res <- func(X = X, Y = Y, ncomp = ncomp, keepX = keepX, design = design, study = HMA.MAE$Batch)
      } else {
        L4 <- c("\nPerforming mint.block.spls")
        cat(L4)
        res <- func(X = X, Y = Y, ncomp = ncomp, keepX = keepX, design = design, keepY = keepY, mode = mode, study = HMA.MAE$Batch)
      }
    } else {
      if (DA){
        L4 <- c("\nPerforming block.splsda")
        cat(L4)
        res <- func(X = X, Y = Y, ncomp = ncomp, keepX = keepX, design = design)
      } else {
        L4 <- c("\nPerforming block.spls")
        cat(L4)
        res <- func(X = X, Y = Y, ncomp = ncomp, keepX = keepX, design = design, keepY = keepY, mode = mode)
      }
      
    }
    L5 <- c("\nFinished model, saving results to: ", fname.out)
    cat(L5)
    # Save results
    saveRDS(res, file = fname.out)
    ## Save info/options used to generate output
    L6 <- c("\nOther params chosen...\nmode:", mode,"\npct_cells_detected:",pct_cells_detected,"\nn_hvrs:",n_hvrs,"\nncomp:",ncomp,"\nkeepXs:", keepXs, "\nkeepY:",keepY)
    fileConn<-file(paste0(tools::file_path_sans_ext(fname.out),"_params.txt"))
    writeLines(c(L1,L2,L3,L4,L5,L6))
    close(fileConn)
  } else {
    file_exists(fname.out)
  }
  
}
