##########################################################################
# Purpose: R Script to calculate average Met/Acc scores for each block and add to cell_metada info
# Output: updated cell_metada RDS file.
#
# Date: 31.Jul.23
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes:
##
# To do: 
## 
# Websites:
##
##########################################################################
#Load packages
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(dplyr, data.table, stringr, SingleCellExperiment, 
               gsubfn, tibble, SummarizedExperiment, MultiAssayExperiment, rlang)

# Input files
cell_metadata_file = "./2_results/scNMT/NMT/Autosomal/cell_metadata.tsv"
NMT_MAE_file = "./2_results/scNMT/NMT/Autosomal/NMT_MAE.rds"

# output files
output.file = "./2_results/scNMT/NMT/Autosomal/cell_metadata_updated.tsv"


# Read input files
cat("\nReading in cell metada file:", cell_metadata_file)
cell_metadata <- fread(cell_metadata_file)
cat("\nReading in NMT MAE file: ", NMT_MAE_file)
NMT_MAE <- readRDS(NMT_MAE_file)

# Calculate
cat("Reading in cell metada file:", cell_metadata_file)
spls_assays <- c("met_Window3k1k", "met_Promoter", "met_CGI", "met_H3K27ac", "met_H3K4me3",
                 "acc_Window3k1k", "acc_Promoter", "acc_CGI", "acc_H3K27ac", "acc_H3K4me3")

df_averages <- NULL
df_averages <- data.frame(row.names = assay(NMT_MAE, spls_assays[2]) %>%
                            as.data.frame() %>% colnames(.))

cat("\nCaclulating averages across:\n", spls_assays)
tmp <- NULL
for(ass in spls_assays){
  tmp <-  
    assay(NMT_MAE, ass) %>%
    as.data.frame() %>%
    colMeans(., na.rm = T) %>%
    as.data.frame() %>%
    dplyr::rename(.,!!ass := 1) %>%
    mutate(!!ass := sprintf(.[[ass]]*100, fmt = '%#.2f'))
  
  df_averages <- cbind(df_averages,tmp)
  
}

# Merge with cell metadata
df_averages <- df_averages %>% tibble::rownames_to_column("cell_id")
cell_metadata_upd <- left_join(cell_metadata, df_averages, by = c("sample" = "cell_id"))

# Save
cat("\nSaving updated cell metadata file into:\n", output.file)
fwrite(cell_metadata_upd, file = output.file, sep ="\t")

