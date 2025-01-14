##########################################################################
# Purpose: R Script to automate correlation of expression and global methylation level in specified treatment group.
# Output: Correlation results stored in a .csv file and displayed in a png.
#
# Date: 25.Jan.22
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
# To do:
##########################################################################


# Function: To correlate gene expression vs methylation level in:
# i) Specific treatment types,
# ii) all groups (to be coded)
##
# setwd("A:/sean_Burnard/1_HPC_backup/Rmarkdown_TEtranscripts_bulk_colony_sc")
# https://csgillespie.github.io/efficientR/performance.html - webguide for efficient coding in R

## Packages 
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
library(SingleCellExperiment)
library(scran)
# For parallelisation
library(foreach) 
library(doSNOW)
library(doParallel)

############### Practice (to delete)
# setwd("A:/sean_Burnard/1_HPC_backup/Rmarkdown_TEtranscripts_bulk_colony_sc")


correlate_expr_methyl <- function(sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds",
                                 sample.sheet = "1_data/colony/HL60/Sample_analysis_sheet.csv",
                                 treatments = "dac",
                                 output.fig.dir = "./figures/colony/HL60/DAC_high_vs_low/",
                                 output.file.dir = "./2_results/colony/HL60/DAC_high_vs_low/",
                                 batch.correction.method = "mnnCorrected"){

  sce_rds_file <- sce.rds.file
  sample_sheet <- sample.sheet
  batch_correction_method <- batch.correction.method
  
  ###### Create figure folder ########
  
  mkdirs <- function(fp) {
    if(!file.exists(fp)) {
      mkdirs(dirname(fp))
      dir.create(fp)
    }
  } 
  
  mkdirs(output.fig.dir)
  mkdirs(output.file.dir)
  ############################################################################################################
  ######## Function to update DESeq2 table results (with annotation file info) ###############################
  
  anno_file = "./scripts/annos/gencode.v30_gene_annotation_table3.txt"
  
  Update_with_hgnc <- function(res_ensembl){
    
    gene_IDs <- fread(anno_file) 
    #  gene_IDs <- fread(anno_file) # Read this in within the above function for this to work
    #res_upd <- as.data.frame(res_ensembl) %>% tibble::rownames_to_column( "gene.TE")
    res_upd <- as.data.frame(res_ensembl)
    
    if( sum(str_detect(res_upd$gene.TE, "ENSG*.*")) >= 1){
      gene_IDs$Geneid <- gsub("\\..*","",gene_IDs$Geneid)
    }
    
    res_upd <- res_upd %>% dplyr::left_join(dplyr::distinct(gene_IDs, Geneid, .keep_all = T), by = c('gene.TE' = 'Geneid')) # Merges tables by ensemble_id_names to gain hgnc symbols.
    res_upd <- res_upd %>% mutate(hgnc.TE = coalesce(GeneSymbol, gene.TE))
    res_upd <- res_upd %>% mutate_if(is.character, list(~na_if(.,""))) # To convert/ ensure any 'emtpy' cells contain 'NA' and not ''. This is important for the next command (coalesce)
    res_upd$Type <- ifelse(grepl("ENSG", res_upd$gene.TE), "gene", "TE") # Smooth way to create new column containing whether transcript is from a gene or TE
    #res_upd2 <- res_upd
    res_upd <- res_upd %>% dplyr::select(ensembl.TE = gene.TE, hgnc.TE, Chromosome, Start , End, Strand, Type, Description = Class, 2:(ncol(res_upd)-9))
    
    return(res_upd)
    
  }
  
  ############################################################################################################
  
  message(paste0("working on: ", sce_rds_file, "..."))
  
  ### Import count data and initial clean and filter ############################################################
  
  data <- readRDS(sce_rds_file) #read.delim(sce.rds.file, check.names = F)
  data_sample_sheet <-fread(sample_sheet)
  
  # TO REMOVE ##### Chromosome X, Y and MT (autosomal only)
  #X_to_remove <- which(rowData(data)$Chromosome == "X")
  #Not_X_chromosome <- which(rowData(data)$Chromosome != "X")
  #data_filt <- data[Not_X_chromosome,]
  
  # TO SELECT Select ##### 'prefered' chromsomes
  Autosomal_Chrs <- c("1", "2", "3", "4", "5", "6","7","8", "9","10","11", "12","13","14", "15","16","17", "18","19","20","21", "22", NA) # Autosomal only (Not 'X', 'Y', or 'MT') Need NA to keep TEs...
  Autosomal_genes <- which(rowData(data)$Chromosome %in% Autosomal_Chrs)
  data_filt <- data[Autosomal_genes,]
  
  # Filter for only DAC cell
  data_filt_treat <- data_filt[, data_filt$Treatment == treatments]
  
  # filter for selected batch correction assay (or uncorrected 'logcounts')
  data_filt_treat_and_assay <- assay(data_filt_treat, batch_correction_method)
  
  
  
  # Setting up for loop parallelisaiton (if multiple cores are available) ################################################
  # https://www.blasbenito.com/post/02_parallelizing_loops_with_r/
  N_CORES <- parallel::detectCores()
  my.cl <- parallel::makeCluster((N_CORES-1), 
                                 type = "PSOCK") 
  
  if(length(my.cl) > 1) message("Running loop correlation loop in parallel on ", length(my.cl), " cores! \nThis step can take >10 mins. \nIt takes even longer if it's processed on a single core.")
  
  # registerDoSNOW(cl)
  doParallel::registerDoParallel(cl = my.cl)
  # Run foreach parallelism (need to check if it will successfully run on a single core if only detecting 1 core available)
  
  Sys.time()->start;
  y <- NULL;
  y <- foreach(gene = row.names(data_filt_treat_and_assay), .combine = 'rbind', .inorder=TRUE) %dopar% {
    a <- cor.test(data_filt_treat_and_assay[gene,], data_filt_treat$meanMeth)
    #print(paste(gene, " est:", a$estimate, " p=value:", a$p.value))
    tmp <- data.frame("gene.TE" = gene, "cor.val.estimate" = a$estimate,  "p.value" = a$p.value)
  };
  time_taken_foreach_parallel <- print(Sys.time()-start)
  parallel::stopCluster(cl = my.cl)
  
  ##################################################################
  # Add multiple testing correction column
  y.updated <- y %>% arrange(p.value)
  y.updated$p.adj <- p.adjust(y.updated$p.value, method = "BH" )
  ## Update with gene info
  y.updated <- Update_with_hgnc(y.updated)
  
  
  
  
  #################################################
  # Save results
  write.csv(y.updated, file = paste0(output.file.dir,"Methylation_Exp_",batch_correction_method,"_correlation.csv"), row.names = F)
  ###########################################################
  
  # Basic plots
  # histogram to show distribution of Cor val
  ggplot(y.updated, aes(x=cor.val.estimate)) +
    geom_histogram(bins = 50, colour="black", fill="white")
  
  ggsave(filename = paste0(output.fig.dir,"Expr_Meth_corr_est_hist_", treatments,".png"),
         height = 5, width = 4)
  
  # significance vs. Cor val
  ggplot(y.updated, aes(x=cor.val.estimate, y = -log10(p.adj))) +
    geom_point(size = 0.1) +
    geom_hline(yintercept = -log10(0.05), 
               linetype="dotted", 
               color = "black", size=1.0)
  
  ggsave(filename = paste0(output.fig.dir,"Expr_Meth_corr_est_vs_pval_", treatments,".png"),
         height = 5, width = 4)
  
  message("Finished Expr vs methylation correlation and plots on data file: ",sce.rds.file)
  message("\nCorrelation results file saved in: ",output.file.dir)
  message("\nBasic plots save in: ",output.fig.dir)
}

###############################################################################################################################



