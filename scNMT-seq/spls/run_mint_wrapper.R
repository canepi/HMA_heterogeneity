##########################################################################
# Purpose: R Script to source and run 'mint_wrapper.R' for use on the HPC.
# Output: RDS file containing the results from the MixOmics analysis
#
# Date: 30.Mar.23
# Version: v.0.2.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes:
## 1) Adapted from Al's original HPC script
# To do:
## 1)
# Websites:
##
##########################################################################
source("~/NOMeSeq.HMA/src/bash_qsub/mint_wrapper.R")

mint_wrapper(fun = 'mint.block.spls', 
             mae.path = '~/NOMeSeq.HMA/output_current/Autosomal/NMT_MAE_Tx_only.rds', 
             out.dir = '~/NOMeSeq.HMA/output_current/Autosomal/finalised_models/pct10_Y100-100_Tx_only/', 
             assay.blocks.to.measure = c("met_CGI","met_H3K27ac","met_H3K4me3","met_Promoter","met_Window3k1k", "acc_CGI","acc_H3K27ac","acc_H3K4me3","acc_Promoter","acc_Window3k1k"), # Default will use all available blocks/assay. 
             se.blocks.assay = NA, # Filters specific assay within all SingleExperiments used as blocks. Defaults to the first assay if multiple are available and none specified.
             rna.experiment = 'rna',
             rna.assay = 'rna', # Specify which assay within the RNA experiment.
             mode = 'canonical', 
             pct_cells_detected = 10,
             n_hvrs = 0, 
             batch = NA, 
             ncomp = 3, 
             keepXs = c(50,50,50), 
             keepY = c(100,100,50), 
             force = TRUE)
