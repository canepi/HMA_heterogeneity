if(!require(pacman)){install.packages("pacman");require(pacman)}

pacman::p_load(here, 
               knitr, 
               SingleCellExperiment, 
               argparser, 
               mixOmics, 
               ggplot2)

formals(ggsave)$dpi <- 600

# prevent x11 error on HPC
options(bitmapType = 'cairo', device = 'png')
