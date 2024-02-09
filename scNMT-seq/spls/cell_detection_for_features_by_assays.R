##########################################################################
# Purpose: R Script to calculate cell detection stats for different assays/contexts and plot!
# Output: plots for CpG and GpC cell detection by feature and associated rds file.
#
# Date: 22.03.23
# Version: v.0.2.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) 
# To do: 
##
# Websites:
##
##########################################################################
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(ggplot2, 
               knitr, 
               SingleCellExperiment,
               ggplot2,
               rlang,
               Seurat,
               gsubfn,
               SummarizedExperiment,
               MultiAssayExperiment)



cell_detection_by_assay <- function(input.sce = "A:/sean_Burnard/1_Projects/HMA_Heterogeneity/2_results/scNMT/NOMeSeq/NOMeSeq_MAE.rds",
                                    output.file = "./2_results/scNMT/NOMeSeq/NOMeSeq_MAE_cell_detection.rds",
                                    output.fig ="./figures/scNMT/NOMeSeq/cell_detection",
                                    batch = NA, #can select batches using c("b320", "b620")
                                    plot.width = 9, plot.height = 6){

  ## Create output dir (if not already present)
  dir.create(dirname(output.fig), recursive = TRUE, showWarnings = FALSE)

  cat("\nreading the input file: \n",input.sce,"\n")
  NOMeSeq.MAE <- readRDS(input.sce)
  
  # Filtering for selected batch (if required)
  if (any(is.na(batch)) == FALSE) { # subset if batch is provided. If NA. Keeps all samples
    NOMeSeq.MAE <- NOMeSeq.MAE[,NOMeSeq.MAE$Batch %in% batch]
  }
  cat("\nKept batches:", unique(NOMeSeq.MAE$Batch))
  
  assay.names <- names(assays(NOMeSeq.MAE)) #[c(1,3,11,13)] # Testing just a small subset first.

  cat("\ncalculating the detection rates:", as.vector(rbind('\n', assay.names)))

  detections_list <- list()
  for (assay.name in assay.names) {
    cat("\ncalculating detection rates for...", assay.name)
    assay.value <- assay(NOMeSeq.MAE, assay.name)
    thresholds = seq(0,1,length.out = 21)
    detections <- rep(NA, 21)
    for (i in 1:21){
      detections[i] <- sum(rowSums(!is.na(assay.value))/ncol(assay.value) > thresholds[i])/nrow(assay.value)*100
    }
    detections_list[[assay.name]] <- data.frame(assay = assay.name, threshold = round(thresholds*100), detection = round(detections))
  }
  
  all_detections <- Reduce('rbind', detections_list)

  # Save cell detection info
  cat("\nFinished calculating detection rates, saving info into: ", output.file)
  saveRDS(all_detections, file = output.file)

  # Plots ##########################################################################################################################################################
  cat("Creating plots for each of the contexts separately for CpG and GpC")
  ## CpG
  CpG_Detection <-
    ggplot(all_detections[grepl('met_', all_detections$assay),], aes(threshold, detection)) + 
    geom_point(size=2) + 
    theme_classic() + 
    labs(y='% of Features Detecting', x = '% of Cells Detected') + 
    theme(axis.text=element_text(size=12), 
          axis.title=element_text(size=12,face="bold")) + facet_wrap(.~assay)
  
  ## Save
  cat("\nSaving: \n", paste0(output.fig,"_CpG.png"))
  ggsave(CpG_Detection,
         filename = paste0(output.fig,"_CpG.png"), 
         width = plot.width, height = plot.height)

  ## GpC
  GpC_Detection <-
    ggplot(all_detections[grepl('acc_', all_detections$assay),], aes(threshold, detection)) + 
    geom_point(size=2) + 
    theme_classic() + 
    labs(y='% of Features Detecting', x = '% of Cells Detected') + 
    theme(axis.text=element_text(size=12), 
          axis.title=element_text(size=12,face="bold")) + facet_wrap(.~assay)

  ## Save
  cat("\nSaving: \n", paste0(output.fig,"_GpC.png"))
  ggsave(GpC_Detection,
         filename = paste0(output.fig,"_GpC.png"), 
         width = plot.width, height = plot.height)

  cat("Done! Review plots and proceed to next step(s)...")
}










