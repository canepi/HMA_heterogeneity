##########################################################################
# Purpose: R Script to visualise correlation circos plot from MixOmics (s)PLS functions and save underlying corr matrix
# Output: Circos plots (svg) or tab delim file.
#
# Date: 13.07.23
# Version: v.0.4.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) Contains two functions: 'spls_Circos()' for plotting and 'spls_Circos_cor_mat()' to obtain underlying similarity matrix from the mixomics::circosPlot() plot
# To do: 
## 1)
# Websites:
## https://mixomics-users.discourse.group/t/centrate-variables-in-circos-plot/745/4 (For adjusting circos labels etc )
###########################################################################
# Load library
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(mixOmics,ggplot2,clusterProfiler,enrichplot, dplyr, MultiAssayExperiment, SummarizedExperiment, MatrixGenerics, matrixStats,
               knitr, SingleCellExperiment, rlang, Seurat, gsubfn, data.table, cowplot, gridExtra, grid)



spls_Circos <- function(input.file = './2_results/scNMT/NMT/pct10_Y200/mint.block.spls_NA_rna_rna_pctCells_10_nhvrs_0_ncomps_3_mode_canonical.rds',
                        meta.data.file = "./2_results/scNMT/NMT/cell_metadata.tsv",
                        output.fig.dir = './figures/scNMT/NMT/pct10_Y200/circos/',
                        output.fig.name = 'circos_r0.5_C1.pdf',
                        store.cor.mat = NULL, # Currently not used. In future - If 'Y', specify output.file.dir, if different from the dir used for the input.file. It will ignore corr.cutoff threshold to speed up calc.
                        output.file.dir = NULL, # Not used currently
                        output.file.cor.mat = NULL, # Not used currently
                        corr.cutoff = 0.6,
                        comps.to.include = 1, # Number vector i.e. c(1,2)
                        blocks.to.include = NULL, # specify blocks in a vector if only a subset is required from the spls model.
                        width = 9,
                        height = 9){
  
  #### For testing
  #input.file = './2_results/scNMT/NMT/pct10_Y200/mint.block.spls_NA_rna_rna_pctCells_10_nhvrs_0_ncomps_3_mode_canonical.rds'
  #meta.data.file = "./2_results/scNMT/NMT/cell_metadata.tsv"
  #output.fig.dir = './figures/scNMT/NMT/pct10_Y200/circos/'
  #output.fig.name = 'circos_r0.5_C1.pdf'
  #store.cor.mat = "Y" # If 'Y', specify output.file.dir, if different from the dir used for the input.file
  #output.file.dir = './2_results/scNMT/NMT/pct10_Y200/circos/'
  #output.file.cor.mat = "circos_mat_r0.6_C1.tsv"
  #corr.cutoff <- 0.6
  #comps.to.include <- 1
  #blocks.to.include = NULL # specify blocks in a vector if only a subset is required from the spls model.
  #width = 9
  #height = 9

  ## ------------------------------------------------------------------------ ##
  # Create output.dir
  #dir.create(output.file.dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(output.fig.dir, showWarnings = FALSE, recursive = TRUE)

  cat("\nloading input file ...\n", input.file)
  res <- readRDS(input.file)

  cat("\nloading meta data file ...\n", meta.data.file)
  cell_metadata <- fread(meta.data.file, colClasses = 'character')
  cell_metadata <- data.frame(cell_metadata, row.names = 'sample')


  # TODO add indY in wrapper
  res$indY = which(names(res$X) == 'Y')
  # TODO after rebase of circosPlot branch this should not be necessary
  var.names <- lapply(res$X, colnames)
  res$X <- mapply(res$X, names(res$X), FUN = function(mat, block) {
    colnames(mat) <- paste0(colnames(mat), '_', block)
    mat
  },SIMPLIFY = FALSE)


  #Run circosPlot
  save.output.name <- paste0(output.fig.dir, output.fig.name)
  cat("\nPerforming circos plotting, and saving file into: \n", save.output.name)

  pdf(file = save.output.name, width=width, height=height)
  mixOmics::circosPlot(res, #Works for block.pls, block.spls, block.plsda and block.splsda models (not mint integrated models yet)
                       group = cell_metadata[rownames(res$variates$Y),]$Treatment,
                       cutoff=corr.cutoff,
                       only.show = 'rna',
                       Y.name = "rna",
                       blocks = blocks.to.include,
                       comp = comps.to.include,
                       size.labels = 0.6, block.labels.adj = 0.3)
  dev.off()
  cat("Finished and saved! (hopefully...)")
}


spls_Circos_cor_mat <- function(input.file = './2_results/scNMT/NMT/pct10_Y200/mint.block.spls_NA_rna_rna_pctCells_10_nhvrs_0_ncomps_3_mode_canonical.rds',
                                meta.data.file = "./2_results/scNMT/NMT/cell_metadata.tsv",
                                output.file.dir = './2_results/scNMT/NMT/pct10_Y200/circos/',
                                output.file.name = "circos_mat_C1.tsv",
                                corr.cutoff = 1, # Currently this is permanently set to 1. Because plotting is automatic but not saved, and lower thresholds with more connections slows down the process. Future versions will apply a pre-filtered table using this threshold if specified!
                                comps.to.include = 1,
                                blocks.to.include = NULL){ # specify blocks in a vector if only a subset is required from the spls model.

  #### Testing
  #input.file = './2_results/scNMT/NMT/pct10_Y200/mint.block.spls_NA_rna_rna_pctCells_10_nhvrs_0_ncomps_3_mode_canonical.rds'
  #meta.data.file = "./2_results/scNMT/NMT/cell_metadata.tsv"
  #output.file.dir = './2_results/scNMT/NMT/pct10_Y200/circos/'
  #output.file.name = "circos_mat_C1.tsv"
  #corr.cutoff <- 1 # Currently this is permanently set to 1. Because plotting is automatic but not saved, and lower thresholds with more connections slows down the process. Future versions will apply a pre-filtered table using this threshold if specified!
  #comps.to.include <- 1
  #blocks.to.include = NULL # specify blocks in a vector if only a subset is required from the spls model.


  ## ------------------------------------------------------------------------ ##
  # Create output.dir
  dir.create(output.file.dir, showWarnings = FALSE, recursive = TRUE)
  #dir.create(output.fig.dir, showWarnings = FALSE, recursive = TRUE)

  cat("\nloading input file ...\n", input.file)
  res <- readRDS(input.file)

  cat("\nloading meta data file ...\n", meta.data.file)
  cell_metadata <- fread(meta.data.file, colClasses = 'character')
  cell_metadata <- data.frame(cell_metadata, row.names = 'sample')


  # TODO add indY in wrapper
  res$indY = which(names(res$X) == 'Y')
  # TODO after rebase of circosPlot branch this should not be necessary
  var.names <- lapply(res$X, colnames)
  res$X <- mapply(res$X, names(res$X), FUN = function(mat, block) {
    colnames(mat) <- paste0(colnames(mat), '_', block)
    mat
  },SIMPLIFY = FALSE)

  cat("\nCalculating correlation matrix (that underlies/from the mixomics circos plot)")
  corMat <- mixOmics::circosPlot(res, #Works for block.pls, block.spls, block.plsda and block.splsda models (not mint integrated models yet)
                                 group = cell_metadata[rownames(res$variates$Y),]$Treatment,
                                 cutoff=1,
                                 only.show = 'rna',
                                 Y.name = "rna",
                                 blocks = blocks.to.include,
                                 comp = comps.to.include)

  corMat <- as.data.frame(corMat)

  save.output.name <- paste0(output.file.dir, output.file.name)
  cat("\nFinished, saving into: \n", save.output.name)
  fwrite(corMat,file = save.output.name,
         quote = F, sep = "\t",
         row.names = T, col.names = T)

}