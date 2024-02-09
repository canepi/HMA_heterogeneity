##########################################################################
# Purpose: R Script to create MAE from CpG and GpC pre-computed rates
# Output: Seurat corrected RDS and figures
#
# Date: 17.07.23
# Version: v.0.3.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) This has been 'hard coded' to work in a particular directory tree format for input.
## 2) Added option to automatically remove any row with all NAs
## 3) Added option to filter/select specific batches
## 4) Corrected line 112 to 'tmp_SummExp' from 'tmp_SummarizedExperiment'.
## 5) Added automated creation of table with rowRanges info (key) from experiments in saved MAE.
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
               HDF5Array,
               SummarizedExperiment,
               MultiAssayExperiment,
               ggplot2,
               batchelor,
               rlang,
               scater,
               gsubfn)


####################################################

create_NOMeSeq_MAE <- function(input.dir = "A:/Papers_in_preparation/HMA_Heterogeneity/Data/scNMT-seq/Normalised",
                               dir.ignore.pattern = "MyWinM5S25", # Will ignore any directory/file that matches the pattern here. For multiple patterns, use '|' between expressions.
                               CpG.assay = "rate",      # 'rate' = unnormalised and norm.rate = 'normalised'
                               GpC.assay = "norm.rate", # 'rate' = unnormalised and norm.rate = 'normalised + batch corrected.
                               min.CpG = 5,             # Threshold of minimal number of total CpGs at loci. Replaces them with 'NA' if under this value.
                               min.GpC = 20,            # Threshold of minimal number of total GpCs at loci. Replaces them with 'NA' if under this value.
                               batch = NA,              # NA will default to all available 'batches'. Otherwise c('batch1', 'batch2') to select specific batches by name
                               remove.rows.with.all.NAs = TRUE, # Automatically removes any row with all NAs (after applying 'min' filtering).
                               chr.to.keep = NULL, #
                               output.file ="./2_results/scNMT/NOMeSeq/NOMeSeq_MAE.rds"){
  
  ## Create output dir (if not already present)
  dir.create(dirname(output.file), recursive = TRUE, showWarnings = FALSE)
  
  # Chosen assay type ('rate' for unnormalised and rate.norm for normalised and batch corrected)
  CpG.assay <- CpG.assay
  GpC.assay <- GpC.assay

  # Set thresholds
  min_CpG <- min.CpG    # minimum number of counts required for CpG data
  min_GpC <- min.GpC   # minimum number of counts required for GpC data

  # Specify top level directory with HDF5 files
  dirs <- list.dirs(input.dir, recursive = T)
  dirs <- dirs[!grepl(dir.ignore.pattern, dirs)] # Specific folder(s) to ignore based on specified pattern.

  # CpG
  CpG.dirs <- dirs[grep(pattern = "CpG/", x=dirs)]
  CpG.context <- gsub(x=CpG.dirs, pattern=paste0(input.dir,"/CpG/cached-CpG-"), replacement = "") %>%
    gsub(pattern = 'HL60.|_SeqMonk', replacement =  "") %>%
    paste0("met_", .)
  

  # TEST - limit to two (small) only
  #CpG.dirs <- CpG.dirs[c(1,3)]
  #CpG.context <- CpG.context[c(1,3)]
  cat("\nDetected the following contexts to read in and combine:", as.vector(rbind('\n', CpG.context)))
  
  context.SumExp_list <- list()
  rowRanges_table <- data.frame()
  #test
  #context = 1
  for(context in 1:length(CpG.context)) {
    # context = 1 # TEst only 
    context.dir <- CpG.dirs[context]
    cat("\nReading in HDF5 file from...", context.dir)
    se = HDF5Array::loadHDF5SummarizedExperiment(dir=context.dir)
    
    if (any(is.na(batch)) == FALSE) { # subset if batch is provided. If NA. Keeps all samples
      se <-   se[,se$Batch %in% batch]
    }
    cat("Kept batches:", unique(se$Batch))
    
    if(FALSE == is.null(chr.to.keep) ){
      cat("\nFiltering chromosomes to keep only: \n", chr.to.keep)
      rows_to_keep <- which(as.data.table(rowRanges(se))$seqnames %in%  chr.to.keep)
      se <- se[rows_to_keep,]
    }
  
    # Extract relevant assay and filter out sites with (total) counts under threshold
    context.assay <- as.matrix(assay(se, CpG.assay))
    context.assay.totals <-  as.matrix(assay(se, "total")) # Using total counts (of CpG/GpC in positions to for thresholds)
    is.na(context.assay)  <- which(context.assay.totals < min_CpG)   # Replace with NAs for positions with a total count under threshold.
    
    # Create a new Summarised Experiment with 'min' filtered assay only.
    
    tmp_SummExp <- SummarizedExperiment(assay = list("CpG.rate.f"=context.assay), colData = colData(se),
                         rowRanges = rowRanges(se))
    
    # Remove rows with all NAs (this automatically updates/filters the rowRanges info)
    if(remove.rows.with.all.NAs == TRUE){
      tmp_SummExp <- tmp_SummExp[rowSums(is.na(assay(tmp_SummExp)[ , 1:ncol(assay(tmp_SummExp))])) != ncol(assay(tmp_SummExp)), ]
    }
    
    # Add updated Assay/SummarisedExperiment to list()
    context.SumExp_list[[CpG.context[context]]] <- tmp_SummExp
    
    # Create rowRanges key
    rowRanges_table_tmp <- 
      rowRanges(context.SumExp_list[[CpG.context[context]]]) %>% 
      as.data.frame() %>% 
      dplyr::rename("ensembl_ID" = "ID", "chr" = "seqnames") %>% 
      mutate(context = CpG.context[context]) %>% 
      rownames_to_column("context_ID")
    
    rowRanges_table <- rbind(rowRanges_table, rowRanges_table_tmp)
    
    cat("\nFinished context: ", CpG.context[context])
  }
  

  # GpC
  GpC.dirs <- dirs[grep(pattern = "GpC/", x=dirs)]
  GpC.context <- gsub(x=GpC.dirs, pattern="A:/Papers_in_preparation/HMA_Heterogeneity/Data/scNMT-seq/Normalised/GpC/cached-GpC-", replacement = "") %>%
    gsub(pattern = 'HL60.|_SeqMonk',replacement =  "") %>%
    paste0("acc_", .)
  #print(GpC.context)

  # TEST - limit to two (small) only
  #GpC.dirs <- GpC.dirs[c(1,3)]
  #GpC.context <- GpC.context[c(1,3)]
  cat("\nDetected the following contexts to read in and combine:", as.vector(rbind('\n', GpC.context)))

  for(context in 1:length(GpC.context)) {

    context.dir <- GpC.dirs[context]
    cat("\nReading in HDF5 file from...", context.dir)
    se = HDF5Array::loadHDF5SummarizedExperiment(dir=context.dir)
    
    if (any(is.na(batch)) == FALSE) { # subset if batch is provided. If NA. Keeps all samples
      se <-   se[,se$Batch %in% batch]
    }
    cat("\nKept batches:", unique(se$Batch))
    
    if(FALSE == is.null(chr.to.keep) ){
      cat("\nFiltering chromosomes to keep only: \n", chr.to.keep)
      rows_to_keep <- which(as.data.table(rowRanges(se))$seqnames %in%  chr.to.keep)
      se <- se[rows_to_keep,]
    }
  
    # Extract relevant assay and filter out sites with (total) counts under threshold
    context.assay <- as.matrix(assay(se, GpC.assay))
    context.assay.totals <-  as.matrix(assay(se, "total")) # Using total counts (of GpC/GpC in positions to for thresholds)
    is.na(context.assay)  <- which(context.assay.totals < min_GpC)   # Replace with NAs for positions with a total count under threshold.
    
    # Create a new Summarised Experiment with 'min' filtered assay only.
    
    tmp_SummExp <- SummarizedExperiment(assay = list("GpC.rate.f"=context.assay), colData = colData(se),
                                        rowRanges = rowRanges(se))
    
    # Remove rows with all NAs (this automatically updates/filters the rowRanges info)
    if(remove.rows.with.all.NAs == TRUE){
      tmp_SummExp <- tmp_SummExp[rowSums(is.na(assay(tmp_SummExp)[ , 1:ncol(assay(tmp_SummExp))])) != ncol(assay(tmp_SummExp)), ]
    }
  
    # Add updated assay to a new summarised experiment (only containing the one updated assay) and add to list
    context.SumExp_list[[GpC.context[context]]] <- tmp_SummExp
    
    # Create and update rowRanges key
    rowRanges_table_tmp <- 
      rowRanges(context.SumExp_list[[GpC.context[context]]]) %>% 
      as.data.frame() %>% 
      dplyr::rename("ensembl_ID" = "ID", "chr" = "seqnames") %>% 
      mutate(context = GpC.context[context]) %>% 
      rownames_to_column("context_ID")
    
    rowRanges_table <- rbind(rowRanges_table, rowRanges_table_tmp)
    
    cat("\nFinished context: ", GpC.context[context])
  }
  
  # Create MAE from list of SummarisedExperiments
  cat("\nFinished reading in contexts for CpG and GpCs, merging them together now into a single MultiAssayExperiment (MAE).")
  MAE.norm <- MultiAssayExperiment(experiments = context.SumExp_list, colData = colData(context.SumExp_list[[1]])[,1:4])

  # Save
  cat("\nSaving MAE... ", output.file, "\nThis may take a while and be a relatively large file with all contexts together. Tried to minimise/remove unnecessary info.")
  saveRDS(MAE.norm, file = output.file)
  
  ## RowRanges table
  filename.rowRanges.key <- paste0(tools::file_path_sans_ext(output.file), "_rowRanges_Key.txt")
  cat("\nSaving MAE rowRanges Key... ", filename.rowRanges.key, "\nThis contains the genomic coordinates and other associated info for CpG/GpC rates in each context.")
  write.table(rowRanges_table, 
              file = filename.rowRanges.key, quote = FALSE, sep = "\t", row.names = F)
  
  cat("\nDone! (Finally...) \nHave a great day. :)")
}










