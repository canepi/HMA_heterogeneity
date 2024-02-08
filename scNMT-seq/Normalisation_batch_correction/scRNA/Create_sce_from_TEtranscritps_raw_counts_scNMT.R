##########################################################################
# Purpose: R Script to create sce from updated TEtranscripts count tables for the scNMT data (modified QC/meta data sheet)
# Output: sce object in an rds file and meta data saved as a .tsv
#
# Date: 25.Oct.22
# Version: v.0.1.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## This will flexibly import a single or multiple count tables (enter a common basename or specific single file).
## Still requires a single QC file, but will only import metadata for samples present in the count matrix,.
## Need to specify which column in the sample sheet contains the sample names.
# To do:
##########################################################################

## Packages
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(dplyr, stringr,tidyr, readr, data.table, SingleCellExperiment)

create_sce_scNMT <-  function(input.counts = "./1_data/scNMT/RNAseq/TE_and_gene_counts_Updated_and_Raw",
                              sample.sheet = "./1_data/scNMT/QC_Summary_combined_NMT_all_batches.csv",
                              sample.sheet.name.column = "Sample",
                              output.directory =  "./2_results/scNMT/RNAseq"){
  
  ## ------------- read the raw count tables and merge
  files.to.merge <- list.files(dirname(input.counts), pattern = basename(input.counts), full.names = TRUE)

  if(length(files.to.merge) > 1) {
      cat('Detected multiple matching count tables. The following files will be merged: \n', files.to.merge) 
    } else {
      cat('Importing the one raw count table: ', files.to.merge)
  }

  scRNA_data <- fread(files.to.merge[1], sep = '\t')
  
  if(length(files.to.merge) > 1) {
    for(additional_files in files.to.merge[2:length(files.to.merge)]){ # This imports and merges the additional count tables.
  
      additional_count_table <- fread(additional_files, sep = '\t') %>% select(1,Chromosome,Start, 14:ncol(.)) # Retain only samples and gene ID for merging to prevent duplicate column info.

      scRNA_data <- left_join(scRNA_data, additional_count_table, by = c("ensembl.TE" =  "ensembl.TE", "Chromosome" = "Chromosome", "Start" = "Start"))
    }
  } else ("Imported a single file\n")

  ## ------------- checks
  cat('Number of rows', ':\n',nrow(scRNA_data)) # Should be 60165
  cat('sum of duplicate rows', ':\n') 
  sum(duplicated(scRNA_data)) # There will be duplicate rows identfied here if merging tables accidentally created more rows...
  #> 0

  cat('sum of duplicate ensembl gene IDs', ':\n')
  sum(duplicated(scRNA_data$ensembl.TE))

  cat('sum of duplicate hgnc symbol IDs', ':\n')
  sum(duplicated(scRNA_data$hgnc.TE))
  #> 1554

  # Identified duplicate ensembl.ID due to pseudogenes present on both X and Y chromosomes
  # This step concatenates "_Y" on the duplicate ensembl gene for only the second instance which is specifically on Y.
  # Furthermore, TEtranscripts only seems to map these reads to the X copy (the first in the gene order) as it probably can't differentiate the loci/chromosome
  cat("Appending '_Y' onto duplicate gene ensembl IDs identified on the Y chromosome")
  scRNA_data <-  scRNA_data %>% 
  mutate(ensembl.TE = ifelse(ensembl.TE %in% scRNA_data$ensembl.TE[(duplicated(scRNA_data$ensembl.TE))] & Chromosome == "Y", 
                             paste0(ensembl.TE, "_",Chromosome), ensembl.TE ))
  
  cat('sum of duplicate ensembl gene IDs after appending "_Y"', ':\n')
  sum(duplicated(scRNA_data$ensembl.TE))

  cat('sum of duplicate hgnc symbol IDs', ':\n')
  sum(duplicated(scRNA_data$hgnc.TE))

  ## Create raw counts table
  raw_counts <- scRNA_data %>% dplyr::select(1,14:ncol(scRNA_data))
  raw_counts <- as.matrix(raw_counts, rownames = 'ensembl.TE')

  #############################################################################################################################################
  # Cell meta data
  cat(paste0('Importing the sample sheet for meta data: ', sample.sheet))
  sample_sheet <- fread(sample.sheet, sep = ',')

  cat('Matching sample IOs in the samples sheet to the count table and adding meta data from the sample sheet', '...\n')
  ## ------------- cell/feature metadata
  ## cell metadata (need to order the same as the matrix)
  samples_in_count_table <- colnames(raw_counts)
  sample.sheet.name.column
  cell_metadata <- sample_sheet %>% dplyr::filter(get({{sample.sheet.name.column}}) %in% samples_in_count_table) # Using specified column name to filter for samples in count matrix.
  cell_order <- data.frame("ID_order" = samples_in_count_table) # Create data frame with ID order of cells (so order of meta data is correct for sce object)

  cell_metadata <- cell_order %>% dplyr::left_join(cell_metadata, by = c("ID_order" = {{sample.sheet.name.column}})) %>% #
    dplyr::rename({{sample.sheet.name.column}} := "ID_order") # Changes first column name back to the original specified containing the sample ID's

  # gene metadata
  gene_metadata <- scRNA_data %>% dplyr::select(ID = ensembl.TE, hgnc.TE, Name, Family, Class, Family_Class, Type, Chromosome, Start, End, Strand, Length, gene.class)
  gene_metadata <- data.frame(gene_metadata, row.names = 'ID')

  ## ensure we have the metadata for all features
  stopifnot(nrow(gene_metadata) ==  nrow(raw_counts))

  # Write cell metadata
  cat("writing cell metadata ...\n")
  suppressWarnings(dir.create(output.directory, recursive = T))
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
