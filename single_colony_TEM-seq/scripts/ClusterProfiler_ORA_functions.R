##########################################################################
# Purpose: R Script to perform KEGG or MSigDB
# Output: Output is csv tables of results.
#
# Date: 27.Jan.23
# Version: v.0.2.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes:
## Condensed into ORA for i) MSigDB, KEGG and GO.
## This can loop all selected clusters across all category options selected and save results tables.
## Thus far, only 'Perform_ORA_GO' has compareCluster() as a built in option for these functions.
# To do: 
## Add compareCluster() option for 'Perform_ORA_MSigDB' and 'Perform_ORA_KEGG' (is it even possible for KEGG?)
##########################################################################

#Load packages
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(pathview, dplyr, data.table, stringr, SingleCellExperiment, clusterProfiler, enrichplot, ggplot2, gsubfn, cowplot, clustree, tibble)
# Wesbite references:
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-comparecluster.html
# http://guangchuangyu.github.io/2016/11/showcategory-parameter-for-visualizing-comparecluster-output/

  #########################################################################################################################################################################################
  ###################################################### MSigDB ###########################################################################################################################
  #########################################################################################################################################################################################
library(msigdbr)
message('\nif required run "install.packages("msigdbr")" \n')

Perform_ORA_MSigDB <- function(Cluster.run.file = "./figures/colony/HL60/HVGs/Zscored/Kmeans_cluster_list_Zscored.txt",
                               cluster.number.s = c(1:8), # Will loop to perform each test on each cluster.
                               sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds",  # For background
                               batch_correction_method = "mnnCorrected",
                               chr.to.keep = c("1", "2", "3", "4", "5", "6","7","8", "9","10","11", "12","13","14", "15","16","17", "18","19","20","21", "22", NA), # NA represents TEs 
                               Keep.genes.TEs = c("gene"),
                               
                               p.adj.method = "fdr", 
                               msigdb.categories.to.test = c("H","C1","C2","C3","C4","C5","C6"), # Will automatically perform loop to perform all chosen categories.
                               p.adj.threshold = 0.05,                  
                               q.val.threshold = 0.4,          
                               plot.key.colour = "p.adjust",
                               output.fig.dir = "./figures/colony/HL60/HVGs/Zscored/",
                               output.file.dir = "./figures/colony/HL60/HVGs/Zscored/",
                               output.extension.name = "",
                               add.plots = FALSE){

## http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
## https://www.gsea-msigdb.org/gsea/msigdb

# Variables (to go into function)
#Cluster.run.file = "./figures/colony/HL60/HVGs/Zscored/Kmeans_cluster_list_Zscored.txt"
#cluster.number = 5
#sce_used ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds" # For background
#batch_correction_method = "mnnCorrected"
#p.adj.method = "fdr"
#p.adj.threshold = 0.05
#q.val.threshold = 0.4
#plot.key.colour = "p.adjust"
#output.fig.dir = "./figures/colony/HL60/HVGs/Zscored/"
#output.file.dir = "./figures/colony/HL60/HVGs/Zscored/"
#output.extension.name = ""
#sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds"
#chr.to.keep = c("1", "2", "3", "4", "5", "6","7","8", "9","10","11", "12","13","14", "15","16","17", "18","19","20","21", "22", NA) # NA represents TEs 
#Keep.genes.TEs = c("gene")
#msigdb.categories.to.test <- c("H","C1","C2","C3","C4","C5","C6")

#cluster.number.s = c(1:8)


# Create dir (if required)
dir.create(output.fig.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output.file.dir, recursive = TRUE, showWarnings = FALSE)

##################################
  # Background
  message("Reading in sce datafile for the background list: '", sce.rds.file, "'\n(filtering for 'genes' only.)")
  data <- readRDS(sce.rds.file)

  ## Filter for genes and/or TEs
  cat("\nFiltering for genes and/or TEs: ", Keep.genes.TEs)
  data_filt <- data[rowData(data)$Type %in% Keep.genes.TEs]
  ## Filter Chromosomes
  cat("\nFiltering for selected chromosomes: \n", chr.to.keep)
  chr_to_keep_index <- which(rowData(data_filt)$Chromosome %in% chr.to.keep)
  data_filt <- data_filt[chr_to_keep_index,]
  cat("\nChromosomes kept: \n", unique(rowData(data_filt)$Chromosome))

  background.list <- row.names(data_filt)

  ## Hit List
  cat("\nExtrating hit list from: \n", Cluster.run.file, "\nAnd running analysis for hits in cluster(s): ",cluster.number.s)
  Cluster.run <-  fread(Cluster.run.file) %>% as.data.frame(.)
  
  datalist_run = list()
  ORA_results_combined_run = NULL
    
  for(cluster.number in cluster.number.s){
    
    message("\nFiltering and perform MSigDB analysis on run: ",cluster.number)
    
    ensembl_hit_list <- Cluster.run %>% 
      filter(cluster == cluster.number) %>% 
      dplyr::select(matches("gene")) %>% 
      .[["gene"]]
    cat("Filtering for genes only (no TEs)")
    ensembl_hit_list <- ensembl_hit_list[grepl("ENSG", ensembl_hit_list)]
    cat("\nNumber of genes (Ensembl IDs) in the hit list: \n", length(ensembl_hit_list))


  # Running msigdbr 


  
  ## Select category(ies) to run
    cat("\nPerforming ORA using msigbr categories: \n",msigdb.categories.to.test)
  
    datalist_category = list()
    ORA_results_combined = NULL
    
    for(category in msigdb.categories.to.test){
    
      msigdb.categories = category
      msigdb_to_test  <- msigdbr(species = "Homo sapiens", category = msigdb.categories) %>% 
        dplyr::select(gs_name, human_ensembl_gene )

  ## Use enricher to perform ORA on msigdbr database category.
      em <- enricher(ensembl_hit_list, 
                     TERM2GENE=msigdb_to_test, 
                     universe = background.list,
                     pAdjustMethod = p.adj.method,
                     pvalueCutoff = p.adj.threshold,
                     qvalueCutoff = q.val.threshold
                     )

  # Convert to data table and update so it's safe for exporting.
      ORA_results <- em@result
      ORA_results$GeneRatio <- gsub("/","//",ORA_results$GeneRatio)
      ORA_results$BgRatio <- gsub("/","//",ORA_results$BgRatio)
      ORA_results$msigb.category = msigdb.categories
      ORA_results$cluster.number = cluster.number
  
      datalist_category[[category]] <- ORA_results
    }
    
    ORA_results_combined = do.call(rbind, datalist_category)
    
    datalist_run[[cluster.number]] <- ORA_results_combined
  }
  
  ORA_results_combined_runs = do.call(rbind, datalist_run)
  
  
  
  cat("\nWriting output csv file to: \n",
      paste0(output.file.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_msigdb_results",output.extension.name,".csv") )
  
  write.csv(ORA_results_combined_runs, 
            file = paste0(output.file.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_msigdb_results",output.extension.name,".csv"),
            row.names = F)
  
  cat("\nFinished. If plots are required, set 'add.plots = TRUE' and Select specific cluster.run (not just all of them) and MSigDB category.")
  cat("\nTo make sense of the MSigDB categories, check out:\nhttps://www.gsea-msigdb.org/gsea/msigdb \nhttp://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html \n")
  
  message("\nFYI - The table still preserves results with values under the specified pvalue... The pvalue cut-off will only affect/limit what's plotted in the figures.\n")
  
  if(add.plots == TRUE ){
    message("Code to automate plots still needs to be completed....")
  }
  
} 
  
  #########################################################################################################################################################################################
  ###################################################### KEGG ###########################################################################################################################
  #########################################################################################################################################################################################

Perform_ORA_KEGG <- function(Cluster.run.file = "./figures/colony/HL60/HVGs/Zscored/Kmeans_cluster_list_Zscored.txt",
                             cluster.number.s = c(2,3,5),
                             sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds", # For background
                             batch_correction_method = "mnnCorrected",
                             p.adj.method = "fdr",
                             p.adj.threshold = 1,
                             q.val.threshold = 1,
                             plot.key.colour = "p.adjust",
                             output.fig.dir = "./figures/colony/HL60/HVGs/Zscored/",
                             output.file.dir = "./figures/colony/HL60/HVGs/Zscored/",
                             output.extension.name = "",
                             chr.to.keep = c("1", "2", "3", "4", "5", "6","7","8", "9","10","11", "12","13","14", "15","16","17", "18","19","20","21", "22", NA), # NA represents TEs 
                             Keep.genes.TEs = c("gene"),
                             add.plots = FALSE){
    
# Variables (to go into function)
#Cluster.run.file = "./figures/colony/HL60/HVGs/Zscored/Kmeans_cluster_list_Zscored.txt"
#cluster.number.s = c(2,3,5)
#sce_used ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds" # For background
#batch_correction_method = "mnnCorrected"
#p.adj.method = "fdr"
#p.adj.threshold = 0.05
#q.val.threshold = 0.4
#plot.key.colour = "p.adjust"
#output.fig.dir = "./figures/colony/HL60/HVGs/Zscored/"
#output.file.dir = "./figures/colony/HL60/HVGs/Zscored/"
#output.extension.name = ""
#sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds"
#chr.to.keep = c("1", "2", "3", "4", "5", "6","7","8", "9","10","11", "12","13","14", "15","16","17", "18","19","20","21", "22", NA) # NA represents TEs 
#Keep.genes.TEs = c("gene")
  
  # Create dir (if required)
  dir.create(output.fig.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output.file.dir, recursive = TRUE, showWarnings = FALSE)
  
##################################
  # Background
  message("Reading in sce datafile for the background list: '", sce.rds.file, "'\n(filtering for 'genes' only.)")
  data <- readRDS(sce.rds.file)
  
  ## Filter for genes and/or TEs
  cat("\nFiltering for genes and/or TEs: ", Keep.genes.TEs)
  data_filt <- data[rowData(data)$Type %in% Keep.genes.TEs]
  ## Filter Chromosomes
  cat("\nFiltering for selected chromosomes: \n", chr.to.keep)
  chr_to_keep_index <- which(rowData(data_filt)$Chromosome %in% chr.to.keep)
  data_filt <- data_filt[chr_to_keep_index,]
  cat("\nChromosomes kept: \n", unique(rowData(data_filt)$Chromosome))
  background_list <- row.names(data_filt)
  

  # Hit List
  cat("\nExtrating hit list from: \n", Cluster.run.file, "\nAnd running analysis for hits in cluster(s): ", cluster.number.s)
  Cluster.run <-  fread(Cluster.run.file) %>% as.data.frame(.) %>% filter(cluster %in% cluster.number.s)
  message("\nFiltering for genes only (no TEs)")
  Cluster.run <- dplyr::filter(Cluster.run, grepl("ENSG",gene))
  ensembl_hit_list <- Cluster.run$gene
  
  # Convert to entrez ID
  message("Now converting ENSEMBL IDs to ENTEZ IDs for KEGG using biomaRt")
  library('biomaRt')
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  
  ## Convert hit list
  ids <- ensembl_hit_list
  ENTREZID_table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", 'entrezgene_id', "description"),values=ids,mart= mart)
  ### Merge with original cluster table matching by ensembl ID
  Cluster.run_updated <- left_join(Cluster.run, ENTREZID_table, by = c("gene" = "ensembl_gene_id"))
  
  
  ## Convert background list
  ids <- background_list
  ENTREZID_table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", 'entrezgene_id', "description"),values=ids,mart= mart)
  ENTREZID_list <- ENTREZID_table$entrezgene_id # Extract only ENTREZ ID as vector
  ENTREZID_background_list <- ENTREZID_list[!is.na(ENTREZID_list)] # Remove NA/missing vavlues
  
  
  datalist_run = list()
  KEGG_results_combined_run = NULL
    
    for(cluster.number in cluster.number.s){
      
      message("\nFiltering and perform KEGG analysis on run: ",cluster.number)
      
      ENTREZID_hit_list <- Cluster.run_updated %>% 
        filter(cluster == cluster.number) %>% 
        dplyr::select(matches("entrezgene_id")) %>% 
        .[["entrezgene_id"]]
      
      cat("\nNumber of genes (Ensembl IDs) in the hit list: \n", length(ENTREZID_hit_list))
      
      mkk <- clusterProfiler::enrichKEGG(gene = ENTREZID_hit_list,
                                          keyType = 'ncbi-geneid',
                                          organism = 'hsa',
                                          pvalueCutoff = p.adj.threshold,
                                          qvalueCutoff = q.val.threshold,
                                          pAdjustMethod = p.adj.method,
                                          
                                          universe = as.character(ENTREZID_background_list))
      
      
      KEGG_results <- mkk@result
      KEGG_results$GeneRatio <- gsub("/","//",KEGG_results$GeneRatio)
      KEGG_results$BgRatio <- gsub("/","//",KEGG_results$BgRatio)
      KEGG_results$cluster.number = cluster.number
      KEGG_results$KEGG_type = "KEGG_ORA"
      
      datalist_run[[cluster.number]] <- KEGG_results
    }
    # Combine KEGG runs into a single table
    KEGG_results_combined_run = do.call(rbind, datalist_run)

  
  
    for(cluster.number in cluster.number.s){
      
      message("\nFiltering and perform KEGG analysis on run: ",cluster.number)
      
      ENTREZID_hit_list <- Cluster.run_updated %>% 
        filter(cluster == cluster.number) %>% 
        dplyr::select(matches("entrezgene_id")) %>% 
        .[["entrezgene_id"]]
      
      cat("\nNumber of genes (Ensembl IDs) in the hit list: \n", length(ENTREZID_hit_list))
      
      mkk <- clusterProfiler::enrichMKEGG(gene = ENTREZID_hit_list,
                                          keyType = 'ncbi-geneid',
                                          organism = 'hsa',
                                          pvalueCutoff = p.adj.threshold,
                                          qvalueCutoff = q.val.threshold,
                                          pAdjustMethod = p.adj.method,
                                          
                                          universe = as.character(ENTREZID_background_list))
      
      
      KEGG_results <- mkk@result
      KEGG_results$GeneRatio <- gsub("/","//",KEGG_results$GeneRatio)
      KEGG_results$BgRatio <- gsub("/","//",KEGG_results$BgRatio)
      KEGG_results$cluster.number = cluster.number
      KEGG_results$KEGG_type = "mKEGG_ORA"
      
      datalist_run[[cluster.number]] <- KEGG_results
    }
    # Combine mKEGG runs into a single table
    mKEGG_results_combined_run = do.call(rbind, datalist_run)
    # Combine KEGG and mKEGG results
    KEGG_ORA_results_combined <- rbind(KEGG_results_combined_run,mKEGG_results_combined_run)
  
  
    # Save table
    cat("\nWriting output csv file to: \n",
        paste0(output.file.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_KEGG_results",output.extension.name,".csv") )
    
    write.csv(KEGG_ORA_results_combined, 
              file = paste0(output.file.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_KEGG_results",output.extension.name,".csv"),
              row.names = F)
    cat("\nFinished. If plots are required, set 'add.plots = TRUE' and Select specific cluster.run (not just all of them) and MSigDB category.")
    
    if(add.plots == TRUE ){
      message("Code to automate plots still needs to be completed....")
    }
}

    ########################################################################################################################################################################
    # How to convert an individual vector of ENSEMBL IDs to ENTREZ IDs
    
    # Convert to entrez ID
    #message("Now converting ENSEMBL IDs to ENTEZ IDs for KEGG using biomaRt")
    #library('biomaRt')
    #ids <- ensembl_hit_list
    #mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    #ENTREZID_hit_table <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol", 'entrezgene_id', "description"),values=ids,mart= mart)
    #ENTREZID_hit_list <- ENTREZID_hit_table$entrezgene_id
    #ENTREZID_hit_list <- ENTREZID_hit_list[!is.na(ENTREZID_hit_list)]
    ########################################################################################################################################################################
    
    # Cholesterol/steroid biogenesis pathway 'hsa00100' (cluster 5)
    # How to view pathway:
    ## Run KEGG analysis
    #kEGG_2 <- enrichKEGG(gene = ENTREZID_hit_list,
    #                     keyType = 'ncbi-geneid',
    #                     organism     = 'hsa',
    #                     pvalueCutoff = 0.05)
    # Select KEGG pathway of interest
    # KEGG_Pathway <- "hsa00100"  # Cholesterol/steroid biogenesis pathway 'hsa00100' (cluster 5)
    # use 'browseKEGG()' using stored KEGG result and specified pathway. This will show up in a html window.
    #browseKEGG(kEGG_2, KEGG_Pathway)
    ## Looking at cluster 5!! 'Transcriptional misregulation in cancer' (upreg regardless of methylation). Mentions AML in pathway differentiation resistance.
    # browseKEGG(kEGG_5, 'hsa05202')
    
    
    # https://www.r-bloggers.com/2015/12/tutorial-rna-seq-differential-expression-pathway-analysis-with-sailfish-deseq2-gage-and-pathview/
    # https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
    
    

  #########################################################################################################################################################################################
  ###################################################### GO ORA - loop per cluster, or compareCluster #####################################################################################
  #########################################################################################################################################################################################

Perform_ORA_GO <- function(Cluster.run.file = "./figures/colony/HL60/HVGs/Zscored/Kmeans_cluster_list_Zscored.txt",
                           cluster.number.s = c(4,1,8,7,5,3,6,2), # This will also determine the order of the groups in the plots
                           sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds",  # For background
                           batch_correction_method = "mnnCorrected",
                           chr.to.keep = c("1", "2", "3", "4", "5", "6","7","8", "9","10","11", "12","13","14", "15","16","17", "18","19","20","21", "22", NA), # NA represents TEs 
                           Keep.genes.TEs = c("gene"),# GO ORA doesn't use TEs.
                           
                           Compare.Cluster = FALSE, # If TRUE, compareCluster() will be run only one category (first one if multiple entered). If FALSE, enrichGO() will be performed and looped for all selected clusters and GO categories.
                           categories.to.test = c("BP","CC", "MF"), # Compare.Cluster=TRUE will by default take the first category if multiple are entered.
                           p.adj.method = "fdr",
                           p.adj.threshold = 0.05,        
                           q.val.threshold = 0.4,  
                           plot.key.colour = "p.adjust",
                           showCategory.number = 5, # For compareCluster(), this is the min number of 'top hits' per clusters to display.
                           output.fig.dir = "./figures/colony/HL60/HVGs/Zscored/ORA/",
                           output.file.dir = "./figures/colony/HL60/HVGs/Zscored/ORA/",
                           output.extension.name = "",
                           add.plots = FALSE,
                           
                           plot.dims.manual = NA,
                           plotx.height.png = 22, # Adjust as required (some defaults provided in the script)
                           plotx.width.png = 12, 
                           plotx.height.svg = 20,
                           plotx.width.svg = 12, 
                           ploty.height.png = 15,
                           ploty.width.png = 22, 
                           ploty.height.svg = 10,
                           ploty.width.svg = 17){ 
  
  ## http://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html
  ## https://www.gsea-msigdb.org/gsea/msigdb
  #Cluster.run.file = "./figures/colony/HL60/HVGs/Zscored/Kmeans_cluster_list_Zscored.txt"
  #cluster.number.s = c(4,1,8,7,5,3,6,2)
  #sce.rds.file ="2_results/colony/HL60/Batch_Correction/sce_qc2_norm_mnnCorrect_scMerged_2C.rds"  # For background
  #batch_correction_method = "mnnCorrected"
  #chr.to.keep = c("1", "2", "3", "4", "5", "6","7","8", "9","10","11", "12","13","14", "15","16","17", "18","19","20","21", "22", NA) # NA represents TEs 
  #Keep.genes.TEs = c("gene")
  
  #p.adj.method = "fdr"
  #categories.to.test = c("BP") # Can perform a single one, or any of 'BP', 'CC', 'MF'
  #p.adj.threshold = 0.05                  
  #q.val.threshold = 0.4        
  
  #plot.key.colour = "p.adjust"
  #output.fig.dir = "./figures/colony/HL60/HVGs/Zscored/ORA/"
  #output.file.dir = "./figures/colony/HL60/HVGs/Zscored/ORA/"
  #output.extension.name = ""
  #add.plots = FALSE
  #Compare.Cluster = TRUE
  #showCategory.number = 3 # Function default is 5. This is the min number of 'top hits' per clusters to display.
  
  
  # Create dir (if required)
  dir.create(output.fig.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output.file.dir, recursive = TRUE, showWarnings = FALSE)
  
  ##################################
  # Background
  message("Reading in sce datafile for the background list: '", sce.rds.file, "'\n(filtering for 'genes' only.)")
  data <- readRDS(sce.rds.file)
  
  ## Filter for genes and/or TEs
  cat("\nFiltering for genes and/or TEs: ", Keep.genes.TEs)
  data_filt <- data[rowData(data)$Type %in% Keep.genes.TEs]
  ## Filter Chromosomes
  cat("\nFiltering for selected chromosomes: \n", chr.to.keep)
  chr_to_keep_index <- which(rowData(data_filt)$Chromosome %in% chr.to.keep)
  data_filt <- data_filt[chr_to_keep_index,]
  cat("\nChromosomes kept: \n", unique(rowData(data_filt)$Chromosome))
  
  background.list <- row.names(data_filt)
  
  ## Hit List
  cat("\nExtrating hit list from: \n", Cluster.run.file, "\nAnd running analysis for hits in cluster(s): ",cluster.number.s)
  Cluster.run <-  fread(Cluster.run.file) %>% as.data.frame(.)
  
  
  ############ IF Compare.Cluster == FALSE will loop enrichGO() for selected clusters and GO categories ###############################################################################
  if(Compare.Cluster == FALSE){
    message("Running as enrichGO as Compare.Cluster was set to 'FALSE'")
  
    datalist_run = list()
    ORA_results_combined_run = NULL
  
    # Load organism database before loop
    organism = "org.Hs.eg.db"
    message("Loading organism databse: '", organism, "'")
    #BiocManager::install(organism, character.only = TRUE)
    library(organism, character.only = TRUE)
  
    # Run loop for GO ontology and cluster(s)
    for(cluster.number in cluster.number.s){
      message("\nFiltering and perform GO ORA analysis on cluster: ",cluster.number)
    
      ensembl_hit_list <- Cluster.run %>% 
        filter(cluster == cluster.number) %>% 
        dplyr::select(matches("gene")) %>% 
        .[["gene"]]
      cat("Filtering for genes only (no TEs)")
      ensembl_hit_list <- ensembl_hit_list[grepl("ENSG", ensembl_hit_list)]
      cat("\nNumber of genes (Ensembl IDs) in the hit list: \n", length(ensembl_hit_list))
    
      # Running GO ORA 
      ## Run GO
      message("\nPerforming Go analysis...")
    
      datalist_category = list()
      ORA_results_combined = NULL
    
      for(category in categories.to.test){
      
        GO.ont.categories = category
      
        ## Use enricher to perform ORA on msigdbr database category.
        go_enrich <- enrichGO(gene = ensembl_hit_list,
                              universe = background.list,
                              OrgDb = organism,
                              keyType = 'ENSEMBL',
                              readable = T,
                              ont = GO.ont.categories,
                              pAdjustMethod = p.adj.method,
                              pvalueCutoff = p.adj.threshold,
                              qvalueCutoff = q.val.threshold) # Default is 0.2 but tester results didn't meet this threshold.
        
        # Convert to data table and update so it's safe for exporting.
        ORA_results <- go_enrich@result
        ORA_results$GeneRatio <- gsub("/","//",ORA_results$GeneRatio)
        ORA_results$BgRatio <- gsub("/","//",ORA_results$BgRatio)
        ORA_results$category = category
        ORA_results$cluster.number = cluster.number
      
        datalist_category[[category]] <- ORA_results
      }
    
      ORA_results_combined = do.call(rbind, datalist_category)
    
      datalist_run[[cluster.number]] <- ORA_results_combined
    }
  
    ORA_results_combined_runs = do.call(rbind, datalist_run)
  
    cat("\nWriting output csv file to: \n",
        paste0(output.file.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_GO_results",output.extension.name,".csv") )
  
    write.csv(ORA_results_combined_runs, 
              file = paste0(output.file.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_GO_results",output.extension.name,".csv"),
              row.names = F)
  
    cat("\nFinished. If plots are required, set 'add.plots = TRUE' and Select specific cluster.run (not just all of them) and MSigDB category.")
  
    message("\nFYI - The table still preserves results with values under the specified pvalue... The pvalue cut-off will only affect/limit what's plotted in the figures.\n")
  
    if(add.plots == TRUE ){
      message("Code to automate plots still needs to be completed....")
    }
    message("Finished running enrichGO() and saved plots and results. Proceed to next steps or take a break.")
  }
  
  ############ IF Compare.Cluster == TRUE compareClusters will run using the first (or only) GO category specifed #####################################################################
  if(Compare.Cluster == TRUE){
    message("\nRunning as compareCluster() as option was set to 'TRUE'")
    
    ## Filter for selected clusters 
    cat("Filtering for genes only and selected clusters: \n", cluster.number.s)
    Cluster.run <- Cluster.run %>% filter(cluster %in% cluster.number.s) %>% dplyr::filter(., grepl("ENSG",gene))
    
    ## Convert to cluster list
    Cluster.list <- split(Cluster.run$gene, Cluster.run$cluster) 
    
    ## Order clusters
    message("\nOrganising cluster order by the order provided (this will ensure plots follow this order): \n")
    cat(cluster.number.s, "\n")
    Order_of_cluster <- cluster.number.s
    Cluster.list <- Cluster.list[Order_of_cluster]
    
    # levels(Cluster.list) <- Order_of_cluster
    message("List below should display the same specified order. If not, there were missing or mismatching group names from the supplied table.")
    str(Cluster.list)
    
    # GO magic
    cat("Performing compareCluster on GO category: ", categories.to.test[1])
    GO_cluster_results <- compareCluster(Cluster.list, 
                                         universe = background.list,
                                         fun='enrichGO', 
                                         OrgDb='org.Hs.eg.db',
                                         keyType = 'ENSEMBL',
                                         readable = T,
                                         ont = categories.to.test[1],
                                         pAdjustMethod = p.adj.method,
                                         pvalueCutoff = p.adj.threshold,
                                         qvalueCutoff = q.val.threshold)
    
    head(GO_cluster_results)
    
    ## Plot and save (images and results table)
    
    message("\nFinished GO analysis. Saving plots and results table.")
    
    p1 <- clusterProfiler::dotplot(GO_cluster_results, showCategory = showCategory.number)+ 
          theme(axis.text.x = element_text(face="bold", size = 15),
                axis.text.y = element_text(face="bold", size = 15),
                axis.title.x = element_text(size = 13, vjust = -2),
                plot.margin = margin(1,1,1,2, "cm")) 
    
    cat("saving plot with clusters on x-axis to :\n",
        paste0(output.fig.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_compareCluster_GO-",categories.to.test[1],"_plotx.(png)(svg)"))
    
    # Modifying plot dimension based on number of categories. ##########
    if(is.na(plot.dims.manual) == TRUE & showCategory.number == 3){
      plotx.height.png = 22
      plotx.width.png = 12 
      plotx.height.svg = 20
      plotx.width.svg = 12 
      ploty.height.png = 15
      ploty.width.png = 22 
      ploty.height.svg = 10
      ploty.width.svg = 17
    } else if(is.na(plot.dims.manual) == TRUE & showCategory.number == 4){
      plotx.height.png = 22
      plotx.width.png = 12 
      plotx.height.svg = 20
      plotx.width.svg = 12 
      ploty.height.png = 17
      ploty.width.png = 30 
      ploty.height.svg = 15
      ploty.width.svg = 26
    } else if(is.na(plot.dims.manual) == TRUE & showCategory.number > 4){
      plotx.height.png = 24
      plotx.width.png = 12 
      plotx.height.svg = 24
      plotx.width.svg = 12 
      ploty.height.png = 20
      ploty.width.png = 35 
      ploty.height.svg = 21
      ploty.width.svg = 32
    }
    ############################################################
    
    ggsave(p1, filename = paste0(output.fig.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_compareCluster_GO-",categories.to.test[1],"_plotx.png"),
           height = plotx.height.png, width = plotx.width.png)   
    ggsave(p1, filename = paste0(output.fig.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_compareCluster_GO-",categories.to.test[1],"_plotx.svg"),
           height = plotx.height.svg, width = plotx.width.svg)   
    
    
    p2 <- clusterProfiler::dotplot(GO_cluster_results, showCategory = showCategory.number) + 
          coord_flip() + 
          theme(axis.text.x = element_text(angle = 60, hjust=1, size = 18, face="bold"),
                axis.text.y = element_text(face="bold", size = 18),
                axis.title.y = element_text(size = 17, vjust = 3),
                plot.margin = margin(1,1,1,2, "cm"),
                legend.key.size = unit(1.1, 'cm')) 
    
    cat("saving plot  with clusters on y-axis to :\n",
        paste0(output.fig.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_compareCluster_GO-",categories.to.test[1],"_ploty.(png)(svg)"))
    
    ggsave(p2, filename = paste0(output.fig.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_compareCluster_GO-",categories.to.test[1],"_ploty.png"),
           height = ploty.height.png, width = ploty.width.png)
    ggsave(p2, filename = paste0(output.fig.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_compareCluster_GO-",categories.to.test[1],"_ploty.svg"),
           height = ploty.height.svg, width = ploty.width.svg)
    
    cat("Writing GO (",categories.to.test[1],") compareCluster() results :\n",
        paste0(output.file.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_compareCluster_GO-",categories.to.test[1],".csv"), sep ="")
    # Extract table
    ORA_results <- GO_cluster_results@compareClusterResult
    ORA_results$GeneRatio <- gsub("/","//",ORA_results$GeneRatio) # Prevents clashes with excel I/O
    ORA_results$BgRatio <- gsub("/","//",ORA_results$BgRatio)     # Prevents clashes with excel I/O
    
    write.csv(ORA_results, 
              file = paste0(output.file.dir, tools::file_path_sans_ext(basename(Cluster.run.file)), "_compareCluster_GO-",categories.to.test[1],".csv"),
              row.names = F)
    
    message("Finished running compareCluster() and saved plots and results. Proceed to next steps or take a break.")
    
  }
  
} 