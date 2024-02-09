##########################################################################
# Purpose: R Script to combine NOMeSeq MAE with RNAseq assay
# Output: MAE
#
# Date: 17.07.23
# Version: v.0.2.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) 
# To do: 
##
# Websites:
##
##########################################################################
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(ggplot2, dplyr, knitr, SingleCellExperiment, rlang, Seurat, gsubfn)

create_NMT_MAE <- function(input.MAE.NOME ="A:/sean_Burnard/1_Projects/HMA_Heterogeneity/2_results/scNMT/NOMeSeq/NOMeSeq_MAE.rds",
                           input.RNA.sce = "A:/sean_Burnard/1_Projects/HMA_Heterogeneity/2_results/scNMT/RNAseq/Batch_Correction/Seurat/Seurat4_NAs_sce.rds",
                           RNA.assay = "integrated", # Which assay to keep. Others will be dropped in the combined
                           chr.to.keep.RNA = NULL, # set chromosomes with a vector. Applied only to input.RNA.sce object.
                           output.file = "./2_results/scNMT/NMT/NMT_MAE.rds"){

  ## Create output dir (if not already present)
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)

  cat("\nLoading the NOMESeq file...\n", input.MAE.NOME)
  NOMeSeq.MAE <- readRDS(input.MAE.NOME)
  cat("\nLoading the SCE file...\n", input.RNA.sce)
  rna <- readRDS(input.RNA.sce)
  
  if(FALSE == is.null(chr.to.keep.RNA)){
    cat("\nFiltering input.RNA.sce to keep genes only on chromosomes: \n",chr.to.keep.RNA)
    
    rna <- rna[rowData(rna)$Chromosome %in% chr.to.keep.RNA,]
  }
  
  # Update RNA sample(col) names to match NOMeseq.MAE (Treatment.Batch.Well). Easier that way round...
  cat("\nUpdating RNAseq colnames to match NOMeSeq (TX.Batch.WELL)\n")
  colnames(NOMeSeq.MAE)
  colnames(rna)
  update_names <- colData(rna) %>% as.data.frame() %>%
    dplyr::rename("RNAseq_ID" ="Sample") %>%
    mutate(sample = paste(Treatment,Batch,Well, sep = "."), .before = "RNAseq_ID")
  row.names(update_names) <- update_names$sample # Update rownames. Otherwise when they will overwrite the colnames when adding this colData back.

  colData(rna) <- DataFrame(update_names) # This updates colnames and the colData() at the same time!
  colnames(rna) 
  colData(rna)

  # Create new RNA sce retaining only the specified assay (RNA.assay)
  rna_to_add <- SingleCellExperiment(assays = list(rna = assay(rna, RNA.assay)),
                       colData = colData(rna),
                       rowData = rowData(rna))
  
  # initialise MAE
  NMT.MAE <- NOMeSeq.MAE
  # add rna
  cat("\nAdding RNA assay to the MAE\n", "Only adding the assay: '",RNA.assay, "' (storing assay internally as 'rna')", sep = "")
  NMT.MAE <- MultiAssayExperiment(ExperimentList(c(assays(NOMeSeq.MAE), list(rna = rna_to_add))))
  colData(NMT.MAE) # Empty colData in MAE

  # add colData to MAE (using rna data)
  colData(NMT.MAE) <- colData(rna)
  colData(NMT.MAE)
  
  # save
  cat("\nWriting to the file: \n", output.file)
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(NMT.MAE, file = output.file)
}
