##########################################################################
# Purpose: R Script to visualise PCAs from 'calc_pca.R' (NIPALS)
# Output: PCA plots (png)
#
# Date: 30.03.23
# Version: v.0.0.1
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) Added batch_colours options (order given is forced in plot legend)
# To do: 
## 1) 
# Websites:
##
##########################################################################
if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(ggplot2, dplyr, knitr, SingleCellExperiment, rlang, Seurat, gsubfn, data.table, cowplot, gridExtra, grid)

viz_pca <- function(input.folder = "./2_results/scNMT/NOMeSeq/pca/Batch_NA_pctCells_25_nhvrs_5000_ncomps_3",
                    meta.data = "./2_results/scNMT/NMT/cell_metadata.tsv",
                    plot.folder = "./figures/scNMT/NOMeSeq/pca/Batch_NA_pctCells_25_nhvrs_5000_ncomps_3",
                    plot.width = 16,
                    plot.height = 16,
                    treatment_colours = c("Unt" = "#00BA38", "AZA" =  "#619CFF", "DAC"  = "#F8766D"),
                    batch_colours = c("b320" =  "#F8766D", "b620"  = "#00BA38", "b322" = "#619CFF") # This also decides the order
                    ){

  pca.path <- input.folder
  #facet.nrow <- 2

  pca.files <- list.files(pca.path, pattern = '.rds', full.names = TRUE)
  df_list <- list()

  for (pca.file in pca.files) 
    {
    assay.name <- tools::file_path_sans_ext(basename(pca.file))
    pca.res <- readRDS(pca.file)
    df <- as.data.frame(pca.res$scores)
    df$assay <- assay.name
    ## add rownames as column as it will have side-effects  when
    ## rbinding
    df$sample <- rownames(df)
    rownames(df)<- NULL
    df_list[[assay.name]] <- df
    }


  all_df <- Reduce('rbind', df_list)
  all_df <- data.table(all_df)
  length(unique(all_df$sample))

  cat("Reading in metadata from:", meta.data)
  cell_metadata <- fread(meta.data, colClasses = 'character')
  length(unique(cell_metadata$sample))
  setkey(cell_metadata, sample)
  setkey(all_df, sample)
  all_df <- merge.data.table(all_df, cell_metadata, all.y = FALSE)
  
  unique(all_df$Treatment)
  unique(all_df$Batch)

  # plots
  ## Set colours
  treatment_col <- treatment_colours
  batch_col <- batch_colours

  ## Set Batch order by factor
  all_df$Batch <- factor(all_df$Batch, levels = names(batch_colours))


  pcs = list(c(1,2), c(2,3), c(1,3))
  pcs = pcs[[1]]


  for (pcs in list(c(1,2), c(2,3), c(1,3)))
    {
    shape = "Batch"
    group = "Treatment"
    colBy = "Treatment"
    cols = treatment_colours
  
    PCs <- paste0('PC', pcs)
    fig.name <- paste0(c('PCs', pcs), collapse = '')
    ## methylation
    views <- c('met', 'acc')
  
    view <- views[1]
    #colBy <- colBys[1]
  
    for (view in views)
      {
      df <- all_df[grepl(view, assay)]
    
      filename <-  sprintf('%s/%s-%s-%s.png', plot.folder,view, fig.name, "combined")
      dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
    
      assay_plots <- list()
      for(assay_name in unique(df$assay)){
        # Filter for single assay to plot 
        df_filtered <- df %>% filter(assay == assay_name)
      
        p1 <- ggplot(df_filtered, aes_string(PCs[1], PCs[2])) + geom_point(aes_string(col = group, shape = shape)) +
          facet_wrap(.~assay, scales = 'free') +
          scale_colour_manual(values = cols) +
          theme_classic() + labs(col = colBy)     
      
        # Marginal density plot of x (top panel)
        xdensity <- ggplot(df_filtered, aes_string(PCs[1], fill=shape)) + 
          geom_density(alpha=.5)  

        xdensity
        # Marginal density plot of y (right panel)
        ydensity <- ggplot(df_filtered, aes_string(PCs[2], fill=shape)) + 
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
    
      ggsave(cow_p, 
           filename = filename, width = plot.width, height = plot.height, bg = "white")
    
    }
  
  }
}
