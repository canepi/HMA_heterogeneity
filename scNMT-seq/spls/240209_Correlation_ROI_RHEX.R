##### A script to summarise Region-Wise Correlations Computed by Carlos Riveros

### Load Packages

library(data.table)
library(ggplot2)
library(ggpubr)

### Specify directories

wd <- "Z:/Papers_in_preparation/HMA_Heterogeneity/Data/scNMT-seq/RegionWise"
od <- "Z:/Papers_in_preparation/HMA_Heterogeneity/Data/scNMT-seq/RegionWise/Correlations_Summary"

### Load correlations data

setwd(wd)

aa = readRDS(file = 'correlations-CpG-RNA-3.rds')
bb = readRDS(file = 'correlations-GpC-RNA-3.rds')
# cc = readRDS(file = 'correlations-GpC-CpG-3.rds')

### Set up

setwd(od)

Txs <- c("[ALL]", "AZA", "DAC") # Specify which cell groups are of interest: "[ALL]", "AZA", "DAC" and/or "Unt".
#obs <- c(22, 9, 9) # Specify the minimum number of cells (nObs) for a correlation to be considered for each of the cell groups above. 
n.cells<-c(222, 91, 93)

#GOIs <- read.delim(file="GOIs.txt")
#GOIs <- unlist(GOIs)

# GOIs from Dani
GOIs <- c("ENSG00000263961")

print(names(aa[]))
annotations <- c(1, 2, 3, 4, 8, 9, 10) # Choose genomic contexts from the printed list
print(names(aa[annotations]))

### Plot correlations

pdf(file="GOIs.pdf", width = 14, height = 10)


for (g in 1:length(GOIs)){

  ROI_CpG <- as.data.frame(matrix("", ncol = 35)) 
  ROI_GpC <- as.data.frame(matrix("", ncol = 35)) 
  
for (i in annotations){

    tbl=aa[[i]][['cor']]
    f.tbl = tbl[tbl$ID.y==GOIs[g] & tbl$group %in% Txs,]
    colnames(ROI_CpG) <- colnames(f.tbl)
    ROI_CpG<-rbind(ROI_CpG, f.tbl)
    
    tbl=bb[[i]][['cor']]
    f.tbl = tbl[tbl$ID.y==GOIs[g] & tbl$group %in% Txs,]
    colnames(ROI_GpC) <- colnames(f.tbl)
    ROI_GpC<-rbind(ROI_GpC, f.tbl)
    
    }
  
  ROI_CpG <- ROI_CpG[-1,c(3,4,5,8,9,10,13,14,15,17,18,19,20)]
  ROI_CpG[,c(2:5,7:10)] <- apply(ROI_CpG[,c(2:5,7:10)], 2, function(x) as.numeric(x))
  mid.x <- (ROI_CpG$end.x - ROI_CpG$start.x)/2 + ROI_CpG$start.x
  t.cells <- c()
  t.cells[which(ROI_CpG$group==Txs[1])] <- n.cells[1]
  t.cells[which(ROI_CpG$group==Txs[2])] <- n.cells[2]
  t.cells[which(ROI_CpG$group==Txs[3])] <- n.cells[3]
  pc.cells <- ROI_CpG$nObs/t.cells
  ROI_CpG<-cbind(ROI_CpG, mid.x, pc.cells)
  
  if (min(ROI_CpG$cor) < 0){
  
  start <- as.numeric(ROI_CpG$start.y[1])
  end <- as.numeric(ROI_CpG$end.y[1])
  gene <- ROI_CpG$hgnc[1]
  
  gg <- ggplot(ROI_CpG, aes(x=mid.x, y=cor, color=group))
  gg <- gg + geom_point( aes(alpha = -log10(p),size = pc.cells))
  gg <- gg + geom_vline(xintercept=start) 
  gg <- gg + geom_vline(xintercept=end)
  gg <- gg + geom_hline(yintercept= 0.3)
  gg <- gg + geom_hline(yintercept= -0.3)
  gg1 <- gg + ggtitle(paste(names(aa[i]), GOIs[g], gene))
  
  ROI_GpC <- ROI_GpC[-1,c(3,4,5,8,9,10,13,14,15,17,18,19,20)]
  ROI_GpC[,c(2:5,7:10)] <- apply(ROI_GpC[,c(2:5,7:10)], 2, function(x) as.numeric(x))
  mid.x <- (ROI_GpC$end.x - ROI_GpC$start.x)/2 + ROI_GpC$start.x
  t.cells <- c()
  t.cells[which(ROI_GpC$group==Txs[1])] <- n.cells[1]
  t.cells[which(ROI_GpC$group==Txs[2])] <- n.cells[2]
  t.cells[which(ROI_GpC$group==Txs[3])] <- n.cells[3]
  pc.cells <- ROI_GpC$nObs/t.cells
  ROI_GpC<-cbind(ROI_GpC, mid.x, pc.cells)
  
  start <- as.numeric(ROI_GpC$start.y[1])
  end <- as.numeric(ROI_GpC$end.y[1])
  gene <- ROI_GpC$hgnc[1]
  
  gg <- ggplot(ROI_GpC, aes(x=mid.x, y=cor, color=group))
  gg <- gg + geom_point( aes(alpha = -log10(p),size = pc.cells))
  gg <- gg + geom_vline(xintercept=start) 
  gg <- gg + geom_vline(xintercept=end)
  gg <- gg + geom_hline(yintercept= -0.3)
  gg <- gg + geom_hline(yintercept=0.3)
  gg2 <- gg + ggtitle(paste(names(bb[i]), GOIs[g], gene))

  page<-ggarrange(gg1,gg2, ncol = 1, nrow = 2)
  print(page)
  
  } else{}
  
}

dev.off()

