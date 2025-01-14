##########################################################################
# Purpose: R Script to perform MNN correction
# Output: MNN corrected sce object saved as an RDS file
#
# Date: 06.Mar.23
# Version: v.0.2.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: This was originally developed as part of a package build running on unix/HPC by Al, which would allow shell style submission.
## 1) added assay.type for variable input/selection - defaults still to 'logcounts'.
# To do: 
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
               scater)


batch_correction_mnn <- function(input.file = "2_results/colony/HL60/Batch_Correction/sce_qc2_PerBatchNormPCA.rds", 
                                 output.file = "2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect.rds",
                                 assay.name = "mnnCorrected", # assay name output
                                 assay.input = "logcounts", # assay name input
                                 nhvgs = 5000,
                                 force = FALSE){


if (force | ! file.exists(output.file))
{
    output.file <- output.file
    assay.name <- assay.name
    PCA.name <- paste0('PCA_', assay.name)
    
    assay.type <- assay.input

    dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)

    # https://github.com/LTLA/batchelor/blob/master/vignettes/correction.Rmd

    cat("\nloading the input data...\n")
    sce <- readRDS(input.file)


    ## --------------------------- Batch Correction --------------------------- ##
    cat("\ncorrecting batch effects using full MNN ...\n")
    set.seed(401)
    corrected <- mnnCorrect(sce, batch=sce$Batch)
    cat("\nAdding the corrected logcounts to the SCE...\n")
    assay(sce, assay.name) <- assay(corrected, 'corrected')

    cat("\nperforming PCA using corrected data ...\n")
    sce <- runPCA(sce, exprs_values = assay.name, name = PCA.name)
    cat("Saving file: ", output.file)
    saveRDS(sce, output.file)
} else
{
  file.exists(output.file)
}

}
