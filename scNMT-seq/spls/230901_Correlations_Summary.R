library(data.table)

### Specify directories

wd <- "Z:/Papers_in_preparation/HMA_Heterogeneity/Data/scNMT-seq/RegionWise"
od <- "Z:/Papers_in_preparation/HMA_Heterogeneity/Data/scNMT-seq/RegionWise/Correlations_Summary"

### Load data

setwd(wd)

aa = readRDS(file = 'correlations-CpG-RNA-3.rds')
bb = readRDS(file = 'correlations-GpC-RNA-3.rds')
cc = readRDS(file = 'correlations-GpC-CpG-3.rds')

### Get sPLS RNA Feature Clusters

RNA_clusters <- read.csv(file="Heatmap_NMT_C1-2_RNA_rowKmeans3_colKmeans4-ClusterIDs_GENES.csv")
RNA_clusters <- RNA_clusters[, c(2,3,13)]

### Summarise correlations

summary <- data.frame(matrix("", ncol = 10))

Tx <- c("[ALL]", "AZA", "DAC")
min.obs <- c(22,10,10)
min.cor <- 0
p.val <- 0.05

for (c in 0:3){

for (i in (1: length(names(aa[])))){
  
  tbl=aa[[i]][['cor']]
  
  for (l in (1:length(Tx))) {
  
    f.tbl = tbl[c(group == Tx[l] & nObs >= min.obs[l] & seqnames != 'MT'),]  
  
    if(c>0){
      genes <- RNA_clusters$feature_ID[which(RNA_clusters$cluster == c)]
      f.tbl = f.tbl[which(f.tbl$ID.y %in% genes),]
    }
    
    line <- c()
    
    line[1] <- names(aa[i])
    line[2] <- length(f.tbl$cor) # total correlations
    line[3] <- length(which(f.tbl$p < p.val)) # significant correlations
    line[4] <- length(which(f.tbl$cor < -min.cor & f.tbl$p < p.val)) # negative correlations
    line[5] <- length(which(f.tbl$cor > min.cor & f.tbl$p < p.val)) # positive correlations
    line[6] <- length(which(f.tbl$p < p.val))/length(f.tbl$cor) 
    line[7] <- length(which(f.tbl$cor < -min.cor & f.tbl$p < p.val))/length(f.tbl$cor)
    line[8] <- length(which(f.tbl$cor > min.cor & f.tbl$p < p.val))/length(f.tbl$cor)
    line[9] <- Tx[l]
    line[10] <- c
    
    summary <- rbind(summary,line)
  
  }
}

for (j in ( 1: length(names(bb[])) ) ){
  
  tbl=bb[[j]][['cor']]
  
  for (l in (1:length(Tx))) {
    
    f.tbl = tbl[c(group == Tx[l] & nObs >= min.obs[l] & seqnames != 'MT'),]  
    
    if(c>0){
      genes <- RNA_clusters$feature_ID[which(RNA_clusters$cluster == c)]
      f.tbl = f.tbl[which(f.tbl$ID.y %in% genes),]
    }
    
    line <- c()
    
    line[1] <- names(bb[j])
    line[2] <- length(f.tbl$cor) # total correlations
    line[3] <- length(which(f.tbl$p < p.val)) # significant correlations
    line[4] <- length(which(f.tbl$cor < -min.cor & f.tbl$p < p.val)) # negative correlations
    line[5] <- length(which(f.tbl$cor > min.cor & f.tbl$p < p.val)) # positive correlations
    line[6] <- length(which(f.tbl$p < p.val))/length(f.tbl$cor) 
    line[7] <- length(which(f.tbl$cor < -min.cor & f.tbl$p < p.val))/length(f.tbl$cor)
    line[8] <- length(which(f.tbl$cor > min.cor & f.tbl$p < p.val))/length(f.tbl$cor)
    line[9] <- Tx[l]
    line[10] <- c
    
    summary <- rbind(summary,line)
    
  }
  
}}

for (k in ( 1: length(names(cc[])) ) ){
  
  tbl=cc[[k]][['cor']]
  
  for (l in (1:length(Tx))) {
    
    f.tbl = tbl[c(group == Tx[l] & nObs >= min.obs[l] & seqnames != 'MT'),]  
    
    
    line <- c()
    
    line[1] <- names(cc[k])
    line[2] <- length(f.tbl$cor) # total correlations
    line[3] <- length(which(f.tbl$p < p.val)) # significant correlations
    line[4] <- length(which(f.tbl$cor < -min.cor & f.tbl$p < p.val)) # negative correlations
    line[5] <- length(which(f.tbl$cor > min.cor & f.tbl$p < p.val)) # positive correlations
    line[6] <- length(which(f.tbl$p < p.val))/length(f.tbl$cor) 
    line[7] <- length(which(f.tbl$cor < -min.cor & f.tbl$p < p.val))/length(f.tbl$cor)
    line[8] <- length(which(f.tbl$cor > min.cor & f.tbl$p < p.val))/length(f.tbl$cor)
    line[9] <- Tx[l]
    line[10] <- "NA"
    
    summary <- rbind(summary,line)
    
  }
  
}
  


colnames(summary) <- c("Correlation Context", "Correlation Count", "Significant Count", 
                       "Negative Count", "Positive Count", "Significant %", "Negative %", "Positive %",
                       "Group","Gene Cluster")

setwd(od)

write.csv(summary, file="230901_Correlation_Count_Summary.csv")


