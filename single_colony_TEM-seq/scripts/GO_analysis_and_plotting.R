##########################################################################
# Purpose: R Script to perform preliminary GO analysis using 'ClusterProfiler'
# Output: GO results in a csv file and several plots showing the results presented in various ways (png and svg formats)
#
# Date: 05.Sep.23
# Version: v.0.3.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes:
## Added 'show.Category' option. Default is still 15. It will change number of category displayed for i) barplot, ii) dotplot and iii) treeplot 3 
## Added plot.key.colour which can be p.adjust, pvalue or q.value
## Added option for selecting multiple figure extension types
# To do: 
##########################################################################


## Packages 
library(dplyr)
library(plotly)
library(data.table)
library(stringr)
library(SingleCellExperiment)
library(scran)
library(clusterProfiler)
library(enrichplot)





###############


Perform_Go_and_plot <- function(ensembl.hit.list = ensembl_hit_list, # Needs to be a vector. If this is saved in another file, please extract the list/column and make it a vector.
                                background.list = background_list, # Needs to be a vector.
                                p.adj.method = "fdr", # Can be "none", "fdr","bh" etc
                                p.adj.threshold = 0.05, # Default 'always' being 0.05
                                q.val.threshold = 0.4,
                                enrichGO.ontology = "BP",
                                plot.key.colour = "p.adjust", # Can be 'p.adjust', 'pvalue' or 'qvalue'
                                show.Category = 15, # It will change number of categories displayed for i) barplot, ii) dotplot and iii) treeplot 3 
                                output.fig.dir = "./figures/colony/HL60/DAC_high_vs_low/",
                                output.file.dir = "./2_results/colony/HL60/DAC_high_vs_low/",
                                output.extension.name = "", # This can be left blank with "", and all plots will remain with their basic names.
                                figure.extension = c("png", "svg")
                                ){
  # Variables to delete
  #ensembl.hit.list = ensembl_hit_list # Needs to be a vector. If this is saved in another file, please extract the list/column and make it a vector.
  #background.list = background_list # Needs to be a vector.
  #p.adj.method = "fdr" # Can be "NULL", "fdr","bh" etc
  #p.adj.threshold = 0.05 # Default 'always' being 0.05
  #output.fig.dir = "./figures/colony/HL60/DAC_high_vs_low/"
  #output.file.dir = "./2_results/colony/HL60/DAC_high_vs_low/"
  
  
  
  ###### Create figure folder ########
  
  mkdirs <- function(fp) {
    if(!file.exists(fp)) {
      mkdirs(dirname(fp))
      dir.create(fp)
    }
  } 
  
  mkdirs(output.fig.dir)
  mkdirs(output.file.dir)
  
  ##################################
  organism = "org.Hs.eg.db"
  message("Loading organism databse: '", organism, "'")
  
  #BiocManager::install(organism, character.only = TRUE)
  library(organism, character.only = TRUE)
  
  
  
  
  ## Run GO
  message("Performing Go analysis...")
  go_enrich <- enrichGO(gene = ensembl.hit.list,
                        universe = background.list,
                        OrgDb = organism,
                        keyType = 'ENSEMBL',
                        readable = T,
                        ont = enrichGO.ontology,
                        pAdjustMethod = p.adj.method,
                        pvalueCutoff = p.adj.threshold,
                        qvalueCutoff = q.val.threshold) # Default is 0.2 but tester results didn't meet this threshold.
  
  
  message("Finished Go analyis, and now saving results into: '", output.file.dir, "'")
  
  # Replace slash to // to stop excel automatically making it a date...
  # Excel has an issue with any single '/' or ':' even with spaces between the two numbers
  go_results <- go_enrich@result
  go_results$GeneRatio <- gsub("/","//",go_results$GeneRatio)
  go_results$BgRatio <- gsub("/","//",go_results$BgRatio)
  
  
  ## Save tables
  write.csv(go_results, 
            file= if(output.extension.name==""){paste0(output.file.dir,"GO_results.csv")}else{
              paste0(output.file.dir,"GO_results_",output.extension.name,".csv")
            },
            row.names = F)
  write.table(go_results, 
              file= if(output.extension.name==""){paste0(output.file.dir,"GO_results.txt")}else{
                paste0(output.file.dir,"GO_results_",output.extension.name,".txt")}, 
              row.names = F, sep = "\t",
              quote = F)
  
  
  
  
  ### PLOTS ###################################
  message("Starting plots and saving into: '", output.fig.dir, "'")
  
  # Barplot
  barplot(go_enrich,
          drop = TRUE,
          showCategory = show.Category,
          title = "GO Biological Pathways",
          font.size = 8,
          color = plot.key.colour)
  
  for(fig.ext in figure.extension){
    ggsave(filename = 
             if(output.extension.name==""){paste0(output.fig.dir,"EnrichGo_barplot.",fig.ext)}else{
               paste0(output.fig.dir,"EnrichGo_barplot_",output.extension.name,".",fig.ext)}, 
           width = 12.5, height = 11)
  }
  
  
  # Dotplot
  
  enrichplot::dotplot(go_enrich, color = plot.key.colour)
  Go_dotplot <- enrichplot::dotplot(go_enrich, showCategory=show.Category, color = plot.key.colour)
  
  for(fig.ext in figure.extension){
    ggsave(filename = 
             if(output.extension.name==""){paste0(output.fig.dir,"EnrichGo_dotplot.",fig.ext)}else{
               paste0(output.fig.dir,"EnrichGo_dotplot_",output.extension.name,".",fig.ext)}, 
           width = 12.5, height = 11)
  }
  #ggsave(plot=Go_dotplot,
  #      if(output.extension.name==""){paste0(output.fig.dir,"EnrichGo_dotplot.svg")}else{
  #       paste0(output.fig.dir,"EnrichGo_dotplot_",output.extension.name,".svg")},
  #    width=7, height=6.5)
  
  
  
  # Tree plot showing hierarichical clustering of enriched terms.
  edox2 <- pairwise_termsim(go_enrich)
  p1 <- treeplot(edox2, offset = 10, color = plot.key.colour) #offset to increase distance of group bar from labels. Needs to be adjust with ggsave dims too!
  p2 <- treeplot(edox2, hclust_method = "average", offset = 10, color = plot.key.colour) # Offest was larger for this one. But will be different for each result
  # aplot::plot_list(p1, p2, tag_levels='A')
  
  for(fig.ext in figure.extension){
    ggsave(plot = p1, 
           filename = if(output.extension.name==""){paste0(output.fig.dir,"EnrichGo_Tree_plot_1.",fig.ext)}else{
             paste0(output.fig.dir,"EnrichGo_Tree_plot_1_",output.extension.name,".",fig.ext)}, 
           width = 14, height = 10)
  }
  for(fig.ext in figure.extension){
    ggsave(plot = p2, 
           filename = if(output.extension.name==""){paste0(output.fig.dir,"EnrichGo_Tree_plot_2.",fig.ext)}else{
             paste0(output.fig.dir,"EnrichGo_Tree_plot_2_",output.extension.name,".",fig.ext)}, 
           width = 15, height = 15) # To avoid text clash from tree branch and group level, adjust 'offese' (in function) and export dims. Currently this is 'trial and error'
  }
  
  
  treeplot(edox2, offset = 7,
           nCluster = 5, # default is 5
           showCategory = show.Category, # Default is 30
           color = plot.key.colour)
  
  for(fig.ext in figure.extension){
    ggsave(filename = if(output.extension.name==""){paste0(output.fig.dir,"EnrichGo_Tree_plot_3.",fig.ext)}else{
      paste0(output.fig.dir,"EnrichGo_Tree_plot_3_",output.extension.name,".",fig.ext)}, width = 14, height = 10)
  }
  
  treeplot(edox2, offset = 7,
           nCluster = 5, # default is 5
           showCategory = 10, # Default is 30
           color = plot.key.colour)
  
  for(fig.ext in figure.extension){
    ggsave(filename = if(output.extension.name==""){paste0(output.fig.dir,"EnrichGo_Tree_plot_4.",fig.ext)}else{
      paste0(output.fig.dir,"EnrichGo_Tree_plot_4_",output.extension.name,".",fig.ext)}, width = 14, height = 10)
  }
  
  message("Finished GO analysis and plots. \nGo results saved: ", output.file.dir, "\nGo images saved: ", output.fig.dir)
  
}

##################################################################################