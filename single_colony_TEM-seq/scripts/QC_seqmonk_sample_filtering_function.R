##########################################################################
# Purpose: R Script to automate filtering of sample IDs by specified QC thresholds and cell type.
# Output: Produces four text files of samples passing QC:
## i) list of sample IDs passing QC in a single column
## ii) QC thresholds used (also in the column heading of file 'i')
## iii) Breakdown of number of samples passing QC by treatment type.
## iv) Breakdown of number of samples passing QC by replicate and treatment type.
#
# Date: 25.August.22
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
# To do:
##########################################################################

# Load Packages

#library(ggplot2)
library(dplyr)
library(data.table)
#library(stringr)


## Function

Sample_QC_filter <- function(QC.file = "./2_results/colony/QC/QC_stats_seqmonk_hisat2_summarised_ALL_colonies.csv",
                             cell = "HL60",
                             overall.alignment.rate.threshold = 70,
                             Percent.Genes.Measured.threshold = 35,
                             Percent.in_exons.threshold = 65,
                             output.folder = "./2_results/colony/QC"
                              ){

  # Create dir if it doesn't exist
  dir.create(output.folder, showWarnings = FALSE)
  
  ### Load in tables
  
  cat(paste0("reading in '",QC.file, "' \nand filtering for ", cell, " cells."))
  
  QC.file.path <- QC.file
  
  # Filtering for separate cell lines from the single combined QC data table
  QC_combined <- fread(QC.file.path) %>% 
    filter(Cell == cell)
  
  #QC_combined_MOLM <- fread(QC.file.path) %>% 
  #  filter(Cell == "MOLM")
  
  #QC_combined_MV <- fread(QC.file.path) %>% 
  #  filter(Cell == "MV")
  
  ## Set thresholds ################################
  
  
  
  overall_alignment_rate_threshold = overall.alignment.rate.threshold
  Percent_Genes_Measured_threshold = Percent.Genes.Measured.threshold
  Percent_in_exons_threshold = Percent.in_exons.threshold
  
  
  Dataset_to_filter <- QC_combined
  
  # Filter and store results
  Dataset_to_filtered <- 
    Dataset_to_filter %>%
    filter(overall_alignment_rate > overall_alignment_rate_threshold & 
             Percent_Genes_Measured > Percent_Genes_Measured_threshold & 
             Percent_in_exons > Percent_in_exons_threshold) 
  
  Dataset_to_filtered_treatment_summary <- 
    Dataset_to_filtered %>%
    group_by(Treatment) %>%
    summarise(total_passing_QC = n()) %>%
    mutate(Failed = 96 - total_passing_QC,
           Perc_Passing = total_passing_QC/ 96 * 100)
  
  Dataset_to_filtered_rep_summary <- 
    Dataset_to_filtered %>%
    group_by(Treatment, replicate) %>%
    summarise(total_passing_QC = n()) %>%
    mutate(Failed = 32 - total_passing_QC,
           Perc_Passing = total_passing_QC/ 32 * 100)
  
  cat(paste0("For ", Dataset_to_filter$Cell[1], " cells, thresholds were set as: ", "\noverall_alignment_rate_threshold > ", overall_alignment_rate_threshold, "\nPercent_Genes_Measured_threshold > ", Percent_Genes_Measured_threshold, "\nPercent_in_exons_threshold > ", Percent_in_exons_threshold))
  Dataset_to_filtered_treatment_summary
  Dataset_to_filtered_rep_summary
  
  
  
  #### Exporting results
  # Sample IDs only
  #Dataset_to_filtered$Sample
  #fwrite(as.data.frame(Dataset_to_filtered$Sample), 
  #       file = paste0(output.folder, "/Samples_passing_QC.txt"), sep ="\t", col.names = F)
  
  # Summary of Treatments passing
  fwrite(Dataset_to_filtered_treatment_summary, 
         file = paste0(output.folder, "/QC_treatment_passing.txt"), sep = "\t")
  
  # Summary of reps (and Treatments) passing
  fwrite(Dataset_to_filtered_rep_summary, 
         file = paste0(output.folder, "/QC_treatment_and_reps_passing.txt"), sep = "\t")
  
  
  # QC thresholds set
  Thresholds_used_summary <- data.frame(
    Threshold.type = c("overall_alignment_rate", "Percent_Genes_Measured", "Percent_in_exons"),
    Above = c(overall_alignment_rate_threshold, Percent_Genes_Measured_threshold, Percent_in_exons_threshold))
  
  fwrite(Thresholds_used_summary, 
         file = paste0(output.folder, "/QC_thresholds_used.txt"), sep = "\t")
  
  
  # Sample IDs passing QC, and thresholds as header
  details <- paste0("overall_alignment_rate-", overall_alignment_rate_threshold,".Percent_Genes_Measured-",Percent_Genes_Measured_threshold,".Percent_in_exons-",Percent_in_exons_threshold)
  Updated_IDs_header <- setNames(as.data.frame(Dataset_to_filtered$Sample), details)
  
  fwrite(Updated_IDs_header, 
         file = paste0(output.folder, "/Samples_passing_QC.txt"), sep ="\t", col.names = T)
  
  #######
  
  message("\nfinished, and exported output files to '", output.folder, "'")
  
}

