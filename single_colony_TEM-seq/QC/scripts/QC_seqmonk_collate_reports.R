#### Basic script to collate QC summary info of hisat2 and seqmonk reports for colony results.


##### Load packages
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(cowplot)


###### Create figure folders ########

mkdirs <- function(fp) {
  if(!file.exists(fp)) {
    mkdirs(dirname(fp))
    dir.create(fp)
  }
} 

mkdirs("./figures/colony/QC")
mkdirs("./2_results/colony/HL60/QC")
mkdirs("./2_results/colony/MOLM/QC")
mkdirs("./2_results/colony/MV/QC")
mkdirs("./2_results/colony/QC")

### MV ##############################################################################################################################

colony <- "MV"
# Files to import
Seqmonk <- fread(paste0("1_data/colony/",colony, "/Seqmonk/Seqmonk_QC.txt"))
Hisat2 <- fread(paste0("1_data/colony/",colony, "/Multiqc/multiqc_data/multiqc_hisat2.txt"))


Seqmonk$DataStore <- gsub(".bam", "", Seqmonk$DataStore) # Remove extension in sample name

QC_combined <-
  Seqmonk %>%
    left_join(Hisat2, by = c("DataStore" = "Sample")) %>%
    dplyr::rename(Sample = DataStore) %>%
    separate(Sample, c("Well", "Cell", "replicate", "Treatment", "SN"), "_", remove =F) %>%
    separate(Well, into = c("Row", "Column"), sep = "(?<=[A-Za-z])(?=[0-9])", remove =F) %>%
  
    mutate(Row = as.numeric(match(Row, LETTERS)), # This is for plotting the plate layout
           Column = as.numeric(Column))

# Treatment info already in file name!
QC_combined$Treatment <- gsub('[0-9]+', '', QC_combined$Treatment) # This is specific for MV sample names, but shouldn't affect the others as the other cell line IDs as they didn't have TxNumber in their names.
QC_combined$Treatment <- gsub('U', 'unt', QC_combined$Treatment)
QC_combined$Treatment <- gsub('A', 'aza', QC_combined$Treatment)
QC_combined$Treatment <- gsub('D', 'dac', QC_combined$Treatment)

# Storing independently for this colony and arrange columns by preferred order.

QC_combined_MV <- QC_combined %>% 
  
  select(Sample,Treatment, Cell, replicate, Well, Row, Column, 
         Percent_in_gene = 'Percent in Gene', Percent_in_exons = 'Percent in exons', Percent_in_rRNA = 'Percent in rRNA',
         Percent_Genes_Measured = 'Percent Genes Measured', Percentage_of_max_data_size = 'Percentage of max data size', 
         Percent_in_MT ='Percent in MT', Percent_on_sense_strand = 'Percent on sense strand', 
         paired_total, paired_aligned_none, paired_aligned_one, paired_aligned_multi, paired_aligned_discord_one, unpaired_total,
         unpaired_aligned_none, unpaired_aligned_one, unpaired_aligned_multi, overall_alignment_rate
  )


### MOLM #################################################################################################


colony <- "MOLM"
# Files to import
Seqmonk <- fread(paste0("1_data/colony/",colony, "/Seqmonk/Seqmonk_QC.txt"))
Hisat2 <- fread(paste0("1_data/colony/",colony, "/Multiqc/multiqc_data/multiqc_hisat2.txt"))
# Need to import sample sheet to integrate treatment info
sample_sheet <- fread(paste0("1_data/colony/",colony, "/Sample_analysis_sheet.csv"))



Seqmonk$DataStore <- gsub(".bam", "", Seqmonk$DataStore)

QC_combined <-
  Seqmonk %>%
  left_join(Hisat2, by = c("DataStore" = "Sample")) %>%
  dplyr::rename(Sample = DataStore) %>%
  separate(Sample, c("Cell", "replicate", "Well", "SN"), "_", remove =F) %>%
  separate(Well, into = c("Row", "Column"), sep = "(?<=[A-Za-z])(?=[0-9])", remove =F) %>%
  
  mutate(Row = as.numeric(match(Row, LETTERS)), # This is for plotting the plate layout
         Column = as.numeric(Column))

# Treatment info from sample sheet
QC_combined <-
 sample_sheet %>% 
    select(Sample_ID, Treatment = Treatment_group) %>%
    dplyr::right_join(QC_combined, by = c("Sample_ID" = "Sample")) %>%
    rename(Sample = Sample_ID) 
           

# Storing independently for this colony

QC_combined_MOLM <- QC_combined %>% 
  
  select(Sample,Treatment, Cell, replicate, Well, Row, Column, 
         Percent_in_gene = 'Percent in Gene', Percent_in_exons = 'Percent in exons', Percent_in_rRNA = 'Percent in rRNA',
         Percent_Genes_Measured = 'Percent Genes Measured', Percentage_of_max_data_size = 'Percentage of max data size', 
         Percent_in_MT ='Percent in MT', Percent_on_sense_strand = 'Percent on sense strand', 
         paired_total, paired_aligned_none, paired_aligned_one, paired_aligned_multi, paired_aligned_discord_one, unpaired_total,
         unpaired_aligned_none, unpaired_aligned_one, unpaired_aligned_multi, overall_alignment_rate
         )




### HL60 ###############################################################################################################


colony <- "HL60"

Seqmonk <- fread(paste0("1_data/colony/",colony, "/Seqmonk/Seqmonk_QC.txt"))
Hisat2 <- fread(paste0("1_data/colony/",colony, "/Multiqc/multiqc_data/multiqc_hisat2.txt"))
Hisat2 <- Hisat2[-1]
# Need to import sample sheet to integrate treatment info
sample_sheet <- fread(paste0("1_data/colony/",colony, "/Sample_analysis_sheet.csv"))

Seqmonk$DataStore <- gsub(".bam", "", Seqmonk$DataStore)



QC_combined <-
  Seqmonk %>%
  left_join(Hisat2, by = c("DataStore" = "Sample")) %>%
  dplyr::rename(Sample = DataStore) %>%
  separate(Sample, c("Cell", "colony", "replicate", "Well", "SN"), "_", remove =F) %>%
  separate(Well, into = c("Row", "Column"), sep = "(?<=[A-Za-z])(?=[0-9])", remove =F) %>%
  
  mutate(Row = as.numeric(match(Row, LETTERS)), # This is for plotting the plate layout
         Column = as.numeric(Column))


# Treatment info from sample sheet
QC_combined <-
  sample_sheet %>% 
  select(Sample_ID, Treatment = Treatment_group) %>%
  dplyr::right_join(QC_combined, by = c("Sample_ID" = "Sample")) %>%
  rename(Sample = Sample_ID)

  

# Storing independently for this colony

QC_combined_HL60 <- QC_combined %>% 
  
  select(Sample,Treatment, Cell, replicate, Well, Row, Column, 
         Percent_in_gene = 'Percent in Gene', Percent_in_exons = 'Percent in exons', Percent_in_rRNA = 'Percent in rRNA',
         Percent_Genes_Measured = 'Percent Genes Measured', Percentage_of_max_data_size = 'Percentage of max data size', 
         Percent_in_MT ='Percent in MT', Percent_on_sense_strand = 'Percent on sense strand', 
         paired_total, paired_aligned_none, paired_aligned_one, paired_aligned_multi, paired_aligned_discord_one, unpaired_total,
         unpaired_aligned_none, unpaired_aligned_one, unpaired_aligned_multi, overall_alignment_rate
  )



# Save QC tables (individually and combined) ####################################################################################
colony <- "HL60"
fwrite(QC_combined_HL60, file = paste0("./2_results/colony/",colony,"/QC/QC_stats_seqmonk_hisat2_summarised.csv"))

colony <- "MOLM"
fwrite(QC_combined_MOLM, file = paste0("./2_results/colony/",colony,"/QC/QC_stats_seqmonk_hisat2_summarised.csv"))

colony <- "MV"
fwrite(QC_combined_MV, file = paste0("./2_results/colony/",colony,"/QC/QC_stats_seqmonk_hisat2_summarised.csv"))

# Combined into a single table

fwrite(rbind(QC_combined_HL60,QC_combined_MOLM,QC_combined_MV),
  file = paste0("./2_results/colony/QC/QC_stats_seqmonk_hisat2_summarised_ALL_colonies.csv"))




# Plotting #########################################################################################################################

## Hisat2 vs paired reads (post trimming)

p1 <- 
  QC_combined_HL60 %>%
    ggplot(aes(x= overall_alignment_rate, y = paired_total, col = Treatment, shape = replicate)) +
    geom_point() + scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100), expand = c(0, 0)) +
    theme_minimal() +
    ggtitle("HL60") +
    theme(plot.title = element_text(hjust = 0.5,face="bold"))

p2 <- 
  QC_combined_MOLM %>%
    ggplot(aes(x= overall_alignment_rate, y = paired_total, col = Treatment, shape = replicate)) +
    geom_point() + scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100), expand = c(0, 0)) +
    theme_minimal() +
    ggtitle("MOLM") +
    theme(plot.title = element_text(hjust = 0.5,face="bold"))

p3 <- 
  QC_combined_MV %>%
    ggplot(aes(x= overall_alignment_rate, y = paired_total, col = Treatment, shape = replicate)) +
    geom_point() + 
    scale_y_continuous(labels = scales::comma) +
    scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100), expand = c(0, 0)) +
    theme_minimal()+
    ggtitle("MV") +
    theme(plot.title = element_text(hjust = 0.5,face="bold"))




# Multi-plot 
Align_vs_reads <-
  plot_grid(p1, p2, p3,
            nrow =3)

# Save
save_plot(plot = Align_vs_reads, "./figures/colony/QC/Align_vs_reads.png", 
          nrow = 1, ncol = 3,
          bg = 'white',
          base_height = 10, base_width = 4)




## MT boxplot
#https://stackoverflow.com/questions/65303141/change-color-of-points-within-group-while-keeping-jitter-using-ggplot


p1 <- 
QC_combined_HL60 %>%
  ggplot(aes(x= interaction(Treatment, replicate), y = Percent_in_MT)) +
  geom_boxplot(aes(col = Treatment)) +
  geom_jitter(aes(col = Treatment)) +
  theme_minimal() +
  ggtitle("HL60") +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  
  ylim(0,50) +
  labs(x = "replicate") +
  facet_grid(.~ replicate, scale = "free_x", switch = "x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(3, "lines"),
        strip.background = element_rect(fill = "transparent"))

p2 <- 
QC_combined_MOLM %>%
  ggplot(aes(x= interaction(Treatment, replicate), y = Percent_in_MT)) +
  geom_boxplot(aes(col = Treatment)) +
  geom_jitter(aes(col = Treatment)) +
  theme_minimal() +
  ggtitle("MOLM") +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  
  ylim(0,50) +
  labs(x = "replicate") +
  facet_grid(.~ replicate, scale = "free_x", switch = "x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(3, "lines"),
        strip.background = element_rect(fill = "transparent"))


p3 <- 
QC_combined_MV %>%
  ggplot(aes(x= interaction(Treatment, replicate), y = Percent_in_MT)) +
  geom_boxplot(aes(col = Treatment)) +
  geom_jitter(aes(col = Treatment)) +
  theme_minimal() +
  ggtitle("MV") +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  
  ylim(0,50) +
  labs(x = "replicate") +
  facet_grid(.~ replicate, scale = "free_x", switch = "x") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing.x = unit(3, "lines"),
        strip.background = element_rect(fill = "transparent"))


# Multi-plot
MT_boxplots <-
  plot_grid(p1, p2, p3,
            nrow =3)

# save
save_plot(plot = MT_boxplots, "./figures/colony/QC/Percent_in_MT.png", 
          nrow = 3, ncol = 1,
          bg = 'white',
          base_width = 10)





## Perc in Exons vs % Genes Measured

p1 <- 
  QC_combined_HL60 %>%
  ggplot(aes(x= Percent_in_exons, y = Percent_Genes_Measured, col = Treatment, shape = replicate)) +
  geom_point() +
  theme_minimal() +
  ggtitle("HL60") +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100), expand = c(0, 0)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100), expand = c(0, 0)) +
  geom_vline(xintercept = 70, linetype="dotted", color = "black") +
  geom_hline(yintercept = 35, linetype="dotted", color = "black") 
    

p2 <- 
  QC_combined_MOLM %>%
  ggplot(aes(x= Percent_in_exons, y = Percent_Genes_Measured, col = Treatment, shape = replicate)) +
  geom_point() +
  theme_minimal() +
  ggtitle("MOLM") +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100), expand = c(0, 0)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100), expand = c(0, 0)) +
  geom_vline(xintercept = 70, linetype="dotted", color = "black") +
  geom_hline(yintercept = 35, linetype="dotted", color = "black") 

p3 <- 
  QC_combined_MV %>%
  ggplot(aes(x= Percent_in_exons, y = Percent_Genes_Measured, col = Treatment, shape = replicate)) +
  geom_point() +
  theme_minimal() +
  ggtitle("MV") +
  theme(plot.title = element_text(hjust = 0.5,face="bold")) +
  scale_x_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100), expand = c(0, 0)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 10), limits = c(0,100), expand = c(0, 0)) +
  geom_vline(xintercept = 70, linetype="dotted", color = "black") +
  geom_hline(yintercept = 35, linetype="dotted", color = "black") 



#Plot_2_save <- 
In_Exons_vs_Genes_Measured <-
  plot_grid(p1, p2, p3,
            nrow =3)

save_plot(plot = In_Exons_vs_Genes_Measured, "./figures/colony/QC/In_Exons_vs_Genes_Measured.png", 
          nrow = 3, ncol = 1,
          bg = 'white')


###################################################################################################################################################