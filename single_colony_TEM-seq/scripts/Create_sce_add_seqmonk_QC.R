##########################################################################
# Purpose: R Script to update sce object and meta data with seqmonk QC results
# Output: i) updated sce object and ii) meta data (.tsv)
#
# Date: 25.Jan.22
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
# To do:
##########################################################################

## Packages
library("dplyr")
library("data.table")
library(SingleCellExperiment)

create_sce_add_seqmonk_qc <-  function(input.sce = "./2_results/colony/HL60/sce.rds",
                                       seqmonk.QC.passing = "./2_results/colony/HL60/QC/Samples_passing_QC.txt",
                                       output.directory =  "./2_results/colony/HL60"){
  
  ## ------------- read in sce object and list of samples passing QC
  
  cat(paste0('\nImporting the sce object: ', input.sce))
  input_sce <- readRDS(input.sce)
  
  cat(paste0('\nImporting list of samples passing QC: ', seqmonk.QC.passing))
  Passing_QC <- fread(seqmonk.QC.passing)
  message('\nNumber of samples detected in list (and passing QC): ', nrow(Passing_QC))
  
  # Identify sce samples that passed QC (using for loop)
  
  # if(colnames(input_sce) %in% as.data.frame(Passing_QC)[,]){"Y"} #Alternative approach
  
  Passed_seqmonk_QC <- NULL
  Passed_seqmonk_QC <- c()
  
  for(name in colnames(input_sce)){
    if(name %in% as.data.frame(Passing_QC)[,])
    {tmp <- "Y"
    Passed_seqmonk_QC <- c(Passed_seqmonk_QC,tmp)} else
    {tmp <- "N"
    Passed_seqmonk_QC <- c(Passed_seqmonk_QC,tmp)}
  }
  
  ## Add cells passing seqmonk QC to sce object
  input_sce$Passed_seqmonk_QC <- Passed_seqmonk_QC
  
  # Extract the updated metadata
  cell_metadata <- as.data.frame(colData(input_sce))
  
  # Update metadata tsv
  ## Check if metadata tsv already exists and print warning
  if(file.exists(paste0(output.directory,'/cell_metadata.tsv')) == "TRUE") message("The file '",output.directory,"/cell_metadata.tsv' already exists.\nTherefore, this file will be replaced with the updated metatadata from '",input.sce,"'." )
  
  fwrite(cell_metadata, file = paste0(output.directory,'/cell_metadata.tsv'))
  
  # Save sce object
  cat("\nUpdating the sce object: '", output.directory,"/sce.rds'", sep ="")
  
  saveRDS(object = input_sce, file = paste0(output.directory,"/sce.rds"))
  
  message('\nFinished updating sce and metadata files in: ',output.directory)
}
