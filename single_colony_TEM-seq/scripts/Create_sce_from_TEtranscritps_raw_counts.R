##########################################################################
# Purpose: R Script to create sce from updated TEtranscripts count tables
# Output: sce object in an rds file and meta data saved as a .tsv
#
# Date: 25.Jan.22
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: This is based on Al's script and modified to work with this dataset and workflow, but aiming to produce the same output for downstream integraiton with his workflow
# To do:
##########################################################################

## Packages
library("dplyr")
library("stringr")
library("tidyr")
library("readr")
library("data.table")
library(SingleCellExperiment)

create_sce_colonies <-  function(input.counts = "./2_results/colony/MOLM_updated/TE_and_gene_counts_Updated_and_Raw.txt",
                                sample.sheet = "./1_data/colony/MOLM/Sample_analysis_sheet_all_reps.csv",
                                output.directory =  "./2_results/colony/MOLM_updated"){

  ## ------------- read the raw counts
  
  cat(paste0('Importing the raw counts table: ', input.counts))
  scRNA_data <- fread(input.counts, sep = '\t')
  
  cat(paste0('Importing the sample sheet for meta data: ', sample.sheet))
  sample_sheet <- fread(sample.sheet, sep = ',')
  
  
  ## ------------- checks
  cat('sum of duplicate rows', ':\n')
  sum(duplicated(scRNA_data))
  #> 0
  
  cat('sum of duplicate gene IDs', ':\n')
  sum(duplicated(scRNA_data$ensembl.TE))
  # Identified duplicate ensembl.ID due to pseudogenes present on both X and Y chromosomes
  # This step concatenates "_Y" on the duplicate ensembl gene for only the second instance which is specifically on Y.
  # Furthermore, TEtranscripts only seems to map these reads to the X copy (the first in the gene order) as it probably can't differentiate the loci/chromosome
  scRNA_data <-  scRNA_data %>% 
    mutate(ensembl.TE = ifelse(ensembl.TE %in% scRNA_data$ensembl.TE[(duplicated(scRNA_data$ensembl.TE))] & Chromosome == "Y", 
                               paste0(ensembl.TE, "_",Chromosome), ensembl.TE ))
  
  #sum(duplicated(scRNA_data2$ensembl.TE)) 
  
  #> 0
  cat('sum of duplicate "Feature"s (should be the same?)', ':\n')
  sum(duplicated(scRNA_data$hgnc.TE))
  #> 1554
  
  
  ## Create raw counts table
  raw_counts <- scRNA_data %>% dplyr::select(1,14:ncol(scRNA_data))
  raw_counts <- as.matrix(raw_counts, rownames = 'ensembl.TE')
  
  
  cat('get cell and feature metadata', '...\n')
  ## ------------- cell/feature metadata
  ## cell metadata (need to order the same as the matrix)
  cell_metadata <- sample_sheet %>% dplyr::select(RNAseq_ID, Sample_ID, Methylation_ID, meanMeth, CpGreads, CpGsites, Batch = Plate, Treatment = Treatment_group, Cells_to_keep = To_be_analysed)
  cell_metadata$RNAseq_ID <- gsub(".bam","", cell_metadata$RNAseq_ID) ## Use RNAseq_ID label (need to remove .bam if present)
  cell_order <- as.data.frame(colnames(raw_counts)) # Obtain ID order of cells (so order of meta data is correct for sce object)
  colnames(cell_order) <- "ID_order"
  
  
  
  
  cell_metadata <- cell_order %>% dplyr::left_join(cell_metadata, by = c("ID_order" = "RNAseq_ID")) %>% 
    dplyr::rename("RNAseq_ID" = "ID_order")  # Bit of a longway round, but better ensures the order is maintained...
  
  ## gene metadata
  gene_metadata <- scRNA_data %>% dplyr::select(ID = ensembl.TE, hgnc.TE, Name, Family, Class, Family_Class, Type, Chromosome, Start, End, Strand, Length, gene.class)
  gene_metadata <- data.frame(gene_metadata, row.names = 'ID')
  
  
  
  ## ensure we have the metadata for all features
  
  stopifnot(nrow(gene_metadata) ==  nrow(raw_counts))
  #gene_metadata <- gene_metadata[rownames(raw_counts),]
  
  #cell_metadata <- data.table(sample = colnames(raw_counts), # lower case 's' to match annos for merge
  #                           Treatment = Treatments,
  #                          Batch = factor(Batches))
  
  
  
  # Write cell metadata
  cat("writing cell metadata ...\n")
  
  file.create(paste0(output.directory,'/cell_metadata.tsv'), showWarnings = FALSE)
  fwrite(cell_metadata, file = paste0(output.directory,'/cell_metadata.tsv'))
  
  # Save sce object
  cat('create SCE object with counts and metadata', '...\n')
  sce <- SingleCellExperiment(assays = list(counts = raw_counts),
                              colData = cell_metadata, #DataFrame(Treatment = Treatments, Batch = Batches),
                              rowData = gene_metadata, #DataFrame(gene_metadata),
                              metadata = paste0("Processed on: ", Sys.Date()))
  
  
  
  cat("writing output file 'sce.rds' ...\n")
  saveRDS(object = sce, file = paste0(output.directory,"/sce.rds"))
}
