##########################################################################
# Purpose: R Script to perform plotVar on spls results
# Output: correlation circle plot (png)
#
# Date: 07.08.23
# Version: v.0.2.0
# Written by: Sean Burnard
# Email: sean.burnard@newcastle.edu.au
# Version notes: 
## 1) 
# To do: 
## 1)
# Websites:
## https://mixomics-users.discourse.group/t/production-of-plots-for-publication/370/4
##########################################################################
# Load packages:

if(!require(pacman)){install.packages("pacman");require(pacman)}
pacman::p_load(mixOmics,ggplot2,clusterProfiler,enrichplot, dplyr, MultiAssayExperiment, SummarizedExperiment, MatrixGenerics, matrixStats,
               knitr, SingleCellExperiment, rlang, Seurat, gsubfn, data.table, cowplot, gridExtra, grid, RColorBrewer)

spls_plotVar <- function(input.spls.file, 
                         comps.to.include = c(1,2), # default will be all. Otherwise specify using a vector 'c(1,2)')
                         output.file.dir, # Where the GO results table will go
                         output.fig.dir,
                         cutoff = 0.5,
                         legend = TRUE,
                         var.names = F,
                         plot.title = NULL, # default is 'Correlation Circle Plot'
                         #figure.extension = "png"
                         ){

  ## ------------------------------------------------------------------------ ##

  cat("loading input file ...\n", input.spls.file)
  res <- readRDS(input.spls.file)
  
  # Create dir
  dir.create(output.fig.dir, showWarnings = FALSE, recursive = TRUE)
  
  ## plotting
  cat("\nPerforming plotVar\n")
  cat("\nAnd saving into:\n", )
  
  blocks.remaining <- res$names$blocks[!res$names$blocks %in% c("Y")]
  blocks.remaining_acc <- blocks.remaining[grepl("acc",blocks.remaining)]
  blocks.remaining_met <- blocks.remaining[grepl("met",blocks.remaining)]
  
  cat("\nFirst with default colours\n")
  
  plotVar(res, 
          comp = comps.to.include,
          var.names = var.names, # Can be a vector. Specifying RNA to show labels only!
          legend = legend, cutoff = cutoff,
          blocks = c("Y", blocks.remaining_met, blocks.remaining_acc),
          # col = c("firebrick", cols.met, cols.acc),
          title = if(TRUE == is.null(plot.title)){"Correlation Circle Plot"} else{plot.title}
          )
  
  # Save
  ggsave(filename =   paste0(output.fig.dir,"spls_plotVar_defaultCol.png"), width = 10, height = 9)
  
  cat("\nNow with modified colours (v1 and v2)\n")
  
  # V1 colours
  pal.met <- colorRampPalette(c("#006CBD", "#3BE2FF"))
  pal.acc <- colorRampPalette(c("#86FF8B", "#006623"))
  cols.met <- pal(length(blocks.remaining_met)) 
  cols.acc <- pal2(length(blocks.remaining_acc))
  
  space_selection <- seq(from = 1, to = (length(blocks.remaining_met)*2)-1,
                         by =(length(blocks.remaining_met)*2)/length(blocks.remaining_met)) # This forces a bigger step between colour selections (by 1)
  cols.met <- pal.met(length(blocks.remaining_met)*2)[space_selection]
  space_selection <- seq(from = 1, to = (length(blocks.remaining_acc)*2)-1,
                         by =(length(blocks.remaining_acc)*2)/length(blocks.remaining_acc))
  cols.acc <- pal.acc(length(blocks.remaining_acc)*2)[space_selection]
  
  
  plotVar(res, 
          comp = comps.to.include,
          var.names = var.names, # Can be a vector. Specifying RNA to show labels only!
          legend = legend, cutoff = cutoff,
          blocks = c("Y", blocks.remaining_met, blocks.remaining_acc),
          col = c("firebrick", cols.met, cols.acc),
          title = if(TRUE == is.null(plot.title)){"Correlation Circle Plot"} else{plot.title}
          )
  
  ggsave(filename =   paste0(output.fig.dir,"spls_plotVar_V1.png"), width = 10, height = 9)
  

  # V2 colours
  pal.met2 <- colorRampPalette(c("#006CBD", "#C15CFF"))
  pal.acc2 <- colorRampPalette(c("#C0DD0B", "#006623"))
  space_selection <- seq(from = 1, to = (length(blocks.remaining_met)*2)-1,
                         by =(length(blocks.remaining_met)*2)/length(blocks.remaining_met)) # This forces a bigger step between colour selections (by 1)
  cols.met <- pal.met2(length(blocks.remaining_met)*2)[space_selection]
  space_selection <- seq(from = 1, to = (length(blocks.remaining_acc)*2)-1,
                         by =(length(blocks.remaining_acc)*2)/length(blocks.remaining_acc))
  cols.acc <- pal.acc2(length(blocks.remaining_acc)*2)[space_selection]
  
  
  plotVar(res, 
          comp = comps.to.include,
          var.names = var.names, # Can be a vector. Specifying RNA to show labels only!
          legend = legend, cutoff = cutoff,
          blocks = c("Y", blocks.remaining_met, blocks.remaining_acc),
          col = c("firebrick", cols.met, cols.acc),
          title = if(TRUE == is.null(plot.title)){"Correlation Circle Plot"} else{plot.title}
          )
  

  ggsave(filename =   paste0(output.fig.dir,"spls_plotVar_V2.png"), width = 10, height = 9)
  


  
}

# Plot var points with names, and varying sizes (instead of circles with dots)
  
#plotVar(res, 
#        var.names = c(TRUE, rep(F, times = (length(blocks.remaining_acc) + length(blocks.remaining_met)))), # Can be a vector. Specifying RNA to show labels only!
#        legend = TRUE, cutoff = 0.5,
#        blocks = c("Y", blocks.remaining_met, blocks.remaining_acc),
#        col = c("brown", cols.met, cols.acc),
#        cex = as.character(c(2,5,4,5,5,5,5,5,5,5,5))
#        )
  
  

# brewer.pal(length(blocks.remaining_met), "BuGn")

#  pal3 <- colorRampPalette(c("#86FF8B", "#006623"))
#  pal3(length(blocks.remaining_acc))
# pal3(length(blocks.remaining_acc))

#  divergent_gradient <- usecol(pal = c("firebrick", "white", "steelblue"), n = 7)[-4]

#  seecol(divergent_gradient, main = "A bi-polar color palette")

# divergent_gradient <- usecol(pal = c("#006CBD", "#88EDFF"), n = 5)
# seecol(divergent_gradient, main = "A bi-polar color palette")

#  cols.met <- divergent_gradient

  
  