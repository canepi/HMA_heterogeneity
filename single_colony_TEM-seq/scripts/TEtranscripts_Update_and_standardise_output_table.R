##########################################################################
# Purpose: Import and update count table from TEtranscripts, using Id conversion table created directly from gencode annotation file used.
# Output: Four text files containing the count tables and additional information for:
## i) genes + TEs 
## and TEs only at the 
## ii) element level (essentially filtering the table for only TEs)
## iii) Family-Class level (sum counts)
## iv) Class level (sum counts)
#
# Date: 25.Mar.21
# Version: v.0.3.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: This has been updated to ensure the output table doesn't contain the ensmemblID.version (only ensemblID with no period).
# This should handle dataframes imported with or without the .version with the ensemblID.
# To do:
##########################################################################

# Steps
# 1)  Import merged count table
# 2)  Update ensemble gene list
# 3)  Summarise family and class counts for the TEs
# 4)  Export updated ensemble gene names for i) original list, ii)


## Packages
library("dplyr")
library("stringr")
library("tidyr")
library("readr")
library("data.table")

## Setwd
# setwd("")
rm(list = ls())

##############################
# Import merged data
#####################



update_table <- function(data, 
                         anno_file = "./scripts/annos/gencode.v30_gene_annotation_table.txt", 
                         output.dir = "./2_results",
                         output.appendix = ""){

  ## Variable data to import
  print(paste0("reading in merged count_table:", data))
  data <-fread(data, header = T, check.names = F)
  ####
  
  data <- rename(data, gene.TE = 'gene/TE')
  head(data)
  
  ## Create column to differentiate between gene and TE
  
  Number_of_Genes <- sum(grepl("ENSG",data$gene.TE)) # Since ensembl IDs are used by TEtranscripts, and no TEs use this nomenclature one can count how many rows start with 'ENSG' to know how many genes are present.
  Number_of_TEs <- (nrow(data) - Number_of_Genes) # The difference between total rows and rows starting with 'ENSG' indicates how many TEs
  
  data <- data %>% mutate(Type = c(rep("gene", Number_of_Genes), rep("TE", Number_of_TEs)) ) # TEtranscripts always outputs genes first then TEs, so it's safe to repeat the number of genes first, then TEs.
  nrow(data)
  
  ### Need to remove 
  # 1) ensembl ID with everything after the period. i.e. ENSG00001923.1 (as the period dictates version or something like that)
  
  #data$gene.TE %>%  gsub("\\..*","",a)
  #data$gene.TE <- str_replace(data$gene.TE,"\\..*","")
  ########## NOW KEEPING fulll ensembl ID (to allow matching with custom ID conversion table made from the annotation files)
  
  # 2) '1_data' and .bam' (in sample IDs) 
  data <- data %>% rename_all(funs(str_replace(., "1_data/", "")))
  data <- data %>% rename_all(funs(str_replace(., ".bam", "")))
  
  ## Sanity check - Should be 60165 rows in total
  ### because TEtranscripts outputs a table with all genes and TEs that it looked for (even if finds '0' counts across all samples)
  nrow(data) == 60165 
  
  ################################## 
  # Updated Gene_ID list!
  #######################################
  ## This conversion file was created directly from the gencode file used for the mapping 'gencode.v30.annotation_NoChr.gtf'
  # This allows for 'perfect' matching between scTE (exports hgnc symbol) and TEtranscripts (exports ensembl ID)
  print(paste0("Using annotation file:", anno_file))
  
  gene_IDs <- fread(anno_file) # read.table(anno_file, header = TRUE, fill = FALSE) OR read_tsv(anno_file, show_col_types = FALSE)
  gene_IDs$GeneSymbol <- trimws(gene_IDs$GeneSymbol, which = c("left"))
  gene_IDs$Geneid_short <-  gsub("\\..*", "", gene_IDs$Geneid) #Create a column to remove 'version' of ensembl ID (period and number)
  
  # Detect if main dataframe being imported has ensemblID.version.
  if (sum(grepl("ENSG*..*", data$gene.TE)) > 1){
    
    dat2 <- data %>% dplyr::left_join(dplyr::distinct(gene_IDs, Geneid, .keep_all = T), by = c('gene.TE' = 'Geneid')) # Merges tables by ensemble_id_names to gain hgnc symbols.
    
    dat2 <- dat2 %>% mutate(hgnc.TE = coalesce(GeneSymbol, gene.TE))
    dat2 <- dat2 %>% mutate(ensembl.TE = coalesce(Geneid_short, gene.TE)) 
    
    dat2 <- dat2 %>% dplyr::select(ensembl.TE, hgnc.TE, Type,  Chromosome, Start, End, Strand, Length, gene.class = Class, 2:(ncol(data)-1))
    
    ## Checking all hgnc and ensembl IDs
    # CHECK2 <- dat2%>%dplyr::select(gene.TE, GeneSymbol, hgnc.TE_IDs )
    
    
  } else if(sum(grepl("ENSG*..*", data$Geneid_short)) == 0) { # No ensembl.ID and version
    
    dat2 <- data %>% dplyr::left_join(dplyr::distinct(gene_IDs, Geneid_short, .keep_all = T), by = c('gene.TE' = 'Geneid_short')) # Merges tables by ensemble_id_names to gain hgnc symbols.
    
    dat2 <- dat2 %>% mutate(hgnc.TE = coalesce(GeneSymbol, gene.TE))
    dat2 <- dat2 %>% mutate(ensembl.TE = coalesce(Geneid_short, gene.TE)) 
    
    dat2 <- dat2 %>% dplyr::select(ensembl.TE = gene.TE, hgnc.TE, Type,  Chromosome, Start, End, Strand, Length, gene.class = Class, 2:(ncol(data)-1))
    
  }
  
  #To_change <- grep("ENSG*..*", data$gene.TE)
  #data$gene.TE[To_change]
  #gsub("\\..*", "", data$gene.TE[To_change])
  
  ##############################################################
  ## Working on TEs only
  #################################################
  
  All_TEs_Original <- filter(dat2, Type == "TE")
  
  All_TEs <- separate(All_TEs_Original, col = ensembl.TE, into = c("Name", "Family", "Class"), sep = ":", remove = F) # Separating TEs by ':' into Name, family and class
  
  All_TEs$Family_Class <- paste(All_TEs$Family, All_TEs$Class, sep = ":")  # creating a column with Family and class combined to enable collapsing same names
  
  
  All_TEs <- All_TEs %>%
    dplyr::select(ensembl.TE, Name, Family, Class, Family_Class, (14:ncol(All_TEs)-1) ) 
  
  ## Extracting sample names
  
  Sample_names <- colnames(All_TEs)[6:ncol(All_TEs)]
  
  ################
  # Summing for all samples!
  
  All_TEs_subset_Class <-
    All_TEs %>% group_by(Class) %>%
    summarise_at(Sample_names, sum)
  
  
  
  All_TEs_subset_Family_Class <-
    All_TEs %>% group_by(Family_Class) %>%
    summarise_at(Sample_names, sum)
  
  
  ### Incorporating TE name, family class info into original table
  
  All_TEs_IDs_Only <- All_TEs %>%dplyr::select(ensembl.TE, Name, Family, Class, Family_Class)
  
  dat3 <- left_join(dat2, All_TEs_IDs_Only, by = c('ensembl.TE'= 'ensembl.TE') )
  
  ## Reorder
  
  dat4 <- dat3 %>% 
    dplyr::select(ensembl.TE, hgnc.TE, Name, Family, Class, Family_Class, (3:(ncol(dat3)-4)) ) 
  
  ###############################################################################################################
  ## Create/ check dirs
  
  mkdirs <- function(fp) {
    if(!file.exists(fp)) {
      mkdirs(dirname(fp))
      dir.create(fp)
    }
  } 
  
  mkdirs(output.dir)
  
  ### Save tables (for later use)
  print(paste0("Saving output files into the directory:", output.dir))
  
  # TEs only
  ## Class level
  write.table(All_TEs_subset_Class, file = paste0(output.dir,"/TE_Class_only_counts", output.appendix,".txt"),
              sep = "\t", quote = F, row.names = F)
  ## Family-Class level
  write.table(All_TEs_subset_Family_Class, file = paste0(output.dir,"/TE_Family_Class_only_counts", output.appendix,".txt"),
              sep = "\t", quote = F, row.names = F)
  ## Element level
  write.table(All_TEs, file = paste0(output.dir,"/TE_only_and_all_counts", output.appendix,".txt"), 
              sep = "\t", quote = F, row.names = F)
  
  # Genes + TEs
  write.table(dat4, file = paste0(output.dir,"/TE_and_gene_counts_Updated_and_Raw", output.appendix,".txt"), 
              sep = "\t", quote = F, row.names = F)
  
  ######## END ####################################################################################################
  


}
