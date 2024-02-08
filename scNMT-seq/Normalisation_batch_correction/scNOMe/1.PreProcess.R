# Pre-process coverage files into sparse coverage matrices
#
# Distilled from analysis on CpG data. This is the 3rd reprocessing of data with the addition of 2022 sequencing batch
# and changed annotation definition and alignment / counting of RNA seq by SB.
#
# CR 20210121, 20220825
#
source('/home/carlos/Projects/BIO19003-201707-Canepi-Lee/docs/20220825-CorrelationAnalysis-3Batches/definitions.R')

setwd(workDir)

library(futile.logger)
library(data.table)
library(Matrix)

source('./methCorUtils.R')

logFile = strftime(Sys.Date(), "preprocess-%Y%m%d.log")
flog.appender(flog.indent_wrapper(appender=appender.tee(logFile), indent="-       "))

.inferColData <- function(cnames, split='_') {
  dd = rbindlist(lapply(strsplit(cnames, split=split), as.list))[, c(2,3,1)]
  colnames(dd) = c('Treatment','Batch','Cell')
  dd[Batch == '060320', Batch := 'b320']
  dd[Batch == '250620', Batch := 'b620']
  dd[Batch == '220331', Batch := 'b322']
  dd$Cell = gsub('^([0-9]+)([A-Z])$', '\\2\\1', dd$Cell)
  dd$Cell = gsub('^([A-Z])([0-9]{1})$', '\\10\\2', dd$Cell)
  return(dd)
}

mvpar = methRegionsPar(logDir=workDir, numWorkers=4, jobname='Pre')

#### Methylation data ######
for(myContext in c('CpG','GpC')) {
  
  flog.info("----\n\t\t   ---- START Preprocess for context %s\n\t\t   ----", myContext)
  
  for(x in names(featureData)) {
    featureData[[x]]$outFile = paste0(workDir, 'methylVariance-', x, '-', myContext, '-unnorm.csv.gz')
  }
  
  covFiles = list.files(path=covDir, pattern=paste0(".*.NOMe.", myContext, ".cov.gz"), full.names=T, recursive=F)

  for(iFD in seq_along(featureData)) {
    FD = featureData[[iFD]]
    featName = names(featureData)[iFD]
    flog.info('>>>>>> Start process of %s %s annotation', myContext, featName)
    
    if(!file.exists(FD$outFile)) {
      
      methRegions = get_regions_file(aFile=FD$dataFile, colInfo=FD$colInfo, zero_based=if(FD$adjustStart) "start" else "none")
      flog.info('Regions: %d', nrow(methRegions))
      ncol.org = ncol(methRegions)
      
      methVar = methVar_dataset(covFiles, methRegions, export_all=T, params=mvpar)
      flog.info('From coverage: %d regions', nrow(methVar))
      cns = colnames(methVar)
      cns = cns[6:length(cns)]
      cnw = sapply(cns, function(t) unlist(strsplit(t,'.',fixed=T))[1], USE.NAMES=F)
      ccd = apply(.inferColData(cnames=cnw), 1, paste0, collapse='.')
      cww = sapply(seq(length(ccd)), function(k) gsub(paste0(cnw[k], '.'), paste0(ccd[k], '_'), cns[k], fixed=T))
      setnames(methVar, old=cns, new=cww)
      
      # This is wrong, need to fix all NA rows
      nn = which(apply(methVar[, ..cww], 1, function(z) all(is.na(z))))
      if(length(nn) > 0) {
        flog.info('Remove rows with all NA:', methVar[nn, .SD, .SDcols=1:5], capture=T)
        methVar = methVar[-nn, ]
      }
      flog.info('Methylation variance data: %d regions', nrow(methVar))
      fwrite(methVar, file=FD$outFile)
      flog.info('Wrote methylation variance to file %s', FD$outFile)
    } else {
      flog.info('Methylation variance file exists, skipped: %s', FD$outFile )
    }
    
  }
}

###### RNASeq Data ######
cacheFile = paste0(workDir, 'rnaseq-data-TE+gene-unnorm.rds')
if(!file.exists(cacheFile)) {
  
  rnaFiles = list.files(path=rnaDir, pattern=".*\\.txt", full.names=T)
  rowData = NULL
  colData = NULL
  countData = NULL
  for(aFile in rnaFiles) {
    seqD = fread(file=aFile)
    flog.info('Processing count file %s:\nnRows:\t%d\nnCols:\t%d\nUniqueID:\t%d\nUniqueID-Chr:\t%d', 
              basename(aFile), nrow(seqD), ncol(seqD), seqD[, uniqueN(ensembl.TE)], uniqueN(seqD[, paste0(ensembl.TE, Chromosome)]))
    nn = grep('[0-9]{,2}[A-H]_[AZADCUntNeg]{3}_[0-9]{6}', colnames(seqD))
    if(length(nn) == 0) {
      nn = grep('[A-H][0-9]{,2}_[AZADCUntNeg]{3}_[0-9]{6}', colnames(seqD))
    }
    seqD[, nid := .N, by='ensembl.TE']
    if(nrow(seqD[nid > 1,]) > 0) {
      flog.info('There are %d duplicate ID entries', nrow(seqD[nid > 1,]))
      # Verify repeated entries
      if(!all(seqD[nid > 1, Chromosome] %in% c('X','Y'))) {
        flog.info('Duplicate ID entries not all on chr X and Y', 
                  seqD[nid > 1 & !Chromosome %in% c('X','Y'), .SD, .SDcols=seq(nn[1]-1)], capture=T)
      }
      if(!all(seqD[nid > 1 & Chromosome == 'Y', ..nn] == 0)) {
        flog.info('Not all duplicate ID entries in chr Y have 0 counts')
      }
      zz = seqD[nid > 1 & Chromosome == 'Y', which=T]
      seqD = seqD[-zz, ]
    }
    aCD = .inferColData(colnames(seqD)[nn])
    # On the assumption that all columns before the first data column are row info....
    aRD = seqD[, .SD, .SDcols=seq(nn[1]-1)]
    if(is.null(rowData)) {
      rowData = copy(aRD)
    } else if(!all.equal(rowData, aRD)) {
      flog.error('Rows data from this RNA seq matrix not the same as previous')
      stop(24)
    }
    colData = rbind(colData, aCD)
    setnames(seqD, old=colnames(seqD)[nn], new=aCD[, paste(Treatment, Batch, Cell, sep='.')])
    countData = cbind(countData, seqD[, ..nn])
  }
  
  # Final checks
  if(uniqueN(colData) != ncol(countData)) {
    flog.error('Information for some columns is repeated')
    stop(25)
  }
  nn = which(apply(countData, 1, function(z) all(z < 1e-5)))
  if(length(nn) > 0) {
    flog.info('Removing %d (out of %d) features with all zero counts, split as:', length(nn), nrow(rowData))
    flog.info('', rowData[nn, .N, by='Type'], capture=T)
    rowData = rowData[-nn, ]
    countData = countData[-nn, ]
  }
  
  # Compute a couple of additional values by row and column
  set(rowData, j=c('AggCount','nSamp.1p','nSamp.5p'), value=list(rowSums(countData, na.rm=T), rowSums(countData > 0), rowSums(countData > 4)))
  nte = rowData[Type == 'TE', which=T]
  set(colData, j=c('nGenes','nTEs'), value=list(colSums(countData[-nte, ] > 0), colSums(countData[nte, ] > 0)))
  rnaseqD = list(rowData=rowData, colData=colData, counts=countData)
  saveRDS(rnaseqD, file=cacheFile)
  flog.info('Saved RNA Seq data to cache file %s', cacheFile)
} else {
  rnaseqD = readRDS(file=cacheFile)
  sapply(rnaseqD, setDT)
  flog.info('Read RNA seq data from cache file %s', cacheFile)
}

stop(0)
