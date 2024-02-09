##########################################################################
# Purpose: R Script to visualise PCAs from MixOmics (s)PLS functions.
# Output: PCA plots (png)
#
# Date: 04.04.2023
# Version: v.0.2.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) 
# To do: 
## 1) 
# Websites:
##
##########################################################################
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(ggplot2, dplyr, knitr, SingleCellExperiment, rlang, Seurat, gsubfn, data.table, cowplot, gridExtra, grid)
source("./scripts/utils_custom.R")


viz_mint_block <- function(input.file = "./2_results/scNMT/NMT/mint.block.spls_NA_rna_rna_pctCells_25_nhvrs_0_ncomps_3_mode_canonical.rds",
                           meta.data = "./2_results/scNMT/NMT/cell_metadata.tsv",
                           output.folder = "./figures/scNMT/NMT",
                           plot.width = 16,
                           plot.height = 16,
                           plot.alpha = 0.8,
                           treatment_colours = c("Unt" = "#00BA38", "AZA" =  "#619CFF", "DAC"  = "#F8766D"),
                           batch_colours = c("b320" =  "#F8766D", "b620"  = "#00BA38", "b322" = "#619CFF") # This also decides the order
                           ){

# Variables to del
#input.file = "./2_results/scNMT/NMT/mint.block.spls_NA_rna_rna_pctCells_25_nhvrs_0_ncomps_3_mode_canonical.rds"
#meta.data = "./2_results/scNMT/NMT/cell_metadata.tsv"
#output.folder = "./figures/scNMT/NMT"
#plot.width = 16
#plot.height = 16
#plot.alpha = 0.8
#treatment_colours = c("Unt" = "#00BA38", "AZA" =  "#619CFF", "DAC"  = "#F8766D")
#batch_colours = c("b320" =  "#F8766D", "b620"  = "#00BA38", "b322" = "#619CFF") # This also decides the order



  # Create output.dir
  dir.create(output.folder, showWarnings = FALSE, recursive = TRUE)



  # Load data
  cat("\nLoading input file: ", input.file)
  res <- readRDS(input.file)
  cat("\nextracting variates ...\n")
  variates <- lapply(res$variates, function(x){
    x <- as.data.table(x, keep.rownames = 'sample')
  })


  if (!('rna' %in% names(variates)))
    {
      ## in spls analysis name Y as 'rna'
      names(variates)[which(names(variates) == 'Y')] <- 'rna'
    }

  variates <- rbind_df_list(variates, new_col = 'assay')
  variates <- as.data.table(variates, keep.rownames = 'sample')
  
  cat("adding cell metadata ...\n")
  cell_metadata <- fread(meta.data, colClasses = 'character')
  variates <- setkey(variates, sample)
  cell_metadata <- setkey(cell_metadata, sample)
  variates <- merge.data.table(variates, cell_metadata, by = 'sample', all.x = TRUE, all.y = FALSE)
  
  comp <- 'comp2'
  # Al previously applied a filter from -10% explained on comp2. Why?
  # variates <- variates[get(comp) > -10]
  assays <- unique(variates$assay)
  if ('rna' %in% assays)
    {
      ## put rna first for ggplot
      ind.rna <- grep("rna", assays)
      assays <- c('rna', assays[-ind.rna])
      variates[,assay := factor(assay, levels = assays, ordered = TRUE)]
    }



  # plots
  ## Set colours
  treatment_col <- treatment_colours
  batch_col <- batch_colours
  
  ## Set Batch order by factor
  variates$Batch <- factor(variates$Batch, levels = names(batch_colours))


  comps = list(c(1,2), c(2,3), c(1,3))
  comps = comps[[1]]
  
  cat("\nPlotting...")
  for (comps in list(c(1,2), c(2,3), c(1,3)))
  {
    cat("\n\nComponents: ", comps, "for...")
    shape = "Batch"
    group = "Treatment"
    colBy = "Treatment"
    cols = treatment_colours
    
    fig.name <- paste0(c('comps', comps), collapse = '')
    comps <- paste0('comp', comps)
    ## methylation
    views <- c('met', 'acc')
    
    view <- views[1]
    #colBy <- colBys[1]
    
    for (view in views)
    {
      cat("\n",view)
      df <- variates[grepl(paste(c('rna',view), collapse = "|"), assay)]      
      
      filename <-  sprintf('%s/%s-%s.png', output.folder,view, fig.name)
      dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
      
      assay_plots <- list()
      for(assay_name in unique(df$assay)){
        # Filter for single assay to plot 
        df_filtered <- df %>% filter(assay == assay_name)
        
        p1 <- ggplot(df_filtered, aes_string(comps[1], comps[2])) + geom_point(aes_string(col = group, shape = shape, alpha = plot.alpha)) +
          facet_wrap(.~assay, scales = 'free') +
          scale_colour_manual(values = cols) +
          theme_classic() + labs(col = colBy)     
        
        # Marginal density plot of x (top panel)
        xdensity <- ggplot(df_filtered, aes_string(comps[1], fill=shape)) + 
          geom_density(alpha=.5)  
        
        xdensity
        # Marginal density plot of y (right panel)
        ydensity <- ggplot(df_filtered, aes_string(comps[2], fill=shape)) + 
          geom_density(alpha=.5, show.legend = FALSE) + coord_flip()
        
        ydensity
        
        blankPlot <- ggplot()+geom_blank(aes(1,1))+
          theme(plot.background = element_blank(), 
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), 
                panel.border = element_blank(),
                panel.background = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_blank(), 
                axis.text.y = element_blank(),
                axis.ticks = element_blank()
          )
        
        Plots_saved = list(scatter = p1, xdensity = xdensity, ydensity = ydensity, blankPlot = blankPlot)
        
        assay_plots[[assay_name]] <- plot_grid(Plots_saved$xdensity, Plots_saved$blankPlot, Plots_saved$scatter, Plots_saved$ydensity, 
                                               ncol=2, nrow=2, rel_widths=c(4, 1.8), rel_heights=c(1.8, 4))
        
      }
      cow_p <- cowplot::plot_grid(plotlist = assay_plots, nrow = 4, ncol = 3)
      
      cat("\nFinished all plots for",view, comps, "saving output to: ", filename)
      ggsave(cow_p, 
             filename = filename, width = plot.width, height = plot.height, bg = "white")
      
    }
    
  }
  cat("\n\nHuzzah!! \nCheck the results (fingers crossed). \nThen... Onto the downstream analysis and continue the fun. ^^")
}