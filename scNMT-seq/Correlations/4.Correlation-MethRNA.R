#
# Perform the methylation correlation analysis
#
# We use the results previously computed and agreed:
# CpG data on 222 cells, rates not batch corrected.
# GpC data on 222 cells, normalised-batch corrected rates
# RNA data on 222 cells, normalised-batch corrected by Sean, end of 2022.
#
# We want to load correlation results per feature and the individual points data corresponding to each pair.
# That means restricting the number of features by the significant part of results.
#
# Data has been computed by scripts and report template notebooks.  I will use the list of feature types
# holding the SummarisedExperiment data with the different assays, but trim it down to only the necessary
# matrices to speed up loading by the shiny app.
#
source('/home/carlos/Projects/BIO19003-201707-Canepi-Lee/docs/20221031-CorrelationAnalysis-3Batches-v2/definitions.R')

setwd(workDir)

library(futile.logger)
library(data.table)
library(Matrix)
library(SummarizedExperiment)

source('./methCorUtils.R')
source('methUtils.R')

logFile = strftime(Sys.Date(), "correlation-%Y%m%d.log")
flog.appender(flog.indent_wrapper(appender=appender.tee(logFile), indent="-       "))
flog.info('==========================\n                   ** Start Correlation Analysis **\n                   ==========================\n')

mvpar = methRegionsPar(logDir=workDir, numWorkers=18, jobname='Corr')

mySampProp = 0.05

mxplDir = paste0(workDir3, 'methXplor/')

BPP = getBPP(params=mvpar)

# Load RNASeq data
rnaFile = paste0(rnaDir, 'Seurat4_NAs_sce.rds')

# Load dataset with normalised and batch corrected (integrated) RNA seq data
# We 
processRNAdataSCE <- function(rnaFile) {
  sce = readRDS(file=rnaFile)
  cd = colData(sce)
  cnames = paste(cd[['Treatment']], cd[['Batch']], cd[['Well']], sep='.')
  colnames(sce) = cnames
  rd = as.data.table(rowData(sce), keep.rownames='ensembl')
  # Discard extra columns
  set(rd, j=c('Name', 'Family', 'Class', 'Family_Class'), value=NULL)
  setnames(rd, old=c('Chromosome', 'Start', 'End', 'Length', 'Strand', 'hgnc.TE', 'gene.class'), 
           new=c('seqnames', 'start', 'end', 'width', 'strand', 'hgnc', 'biotype'))
  ite = which(rd[['Type']] != 'gene')
  rd = rd[-ite, ]
  rd[['idx']] = seq(nrow(rd))
  rnaNormData <<- assay(sce, 'integrated_NAs')[-ite, ]
  rnaCellInfo <<- colData(sce)
  rnaGeneInfo <<- rd
  flog.info('Loaded RNA SCE data from %s, %d cells, %d genes', rnaFile, ncol(rnaNormData), nrow(rnaGeneInfo))
}

processRNAdataSCE(rnaFile)


#############
# CpG - RNA correlation
flog.info('%%\n%%       Process %s Context\n%%', myContext)

myContext = 'CpG'
saveFile = paste0('correlations-', myContext, '-RNA-3d.rds')

myMethAssay = ifelse(myContext == 'CpG', 'rate', 'norm.rate')
myWinTotalRC = ifelse(myContext == 'CpG', 5, 20)
myAnnot = 'Promoter'
myAnnot = 'CGI_SeqMonk'
myAnnot = Sys.getenv('CORR_ANNOT', unset=NA)
if (!is.na(myAnnot)) {
  featureData = featureData[myAnnot]
}

corResults = list()
for (myAnnot in names(featureData)) {
  
  flog.info('\nAnnotation: %s\n', myAnnot)
  dname = paste0(mxplDir, 'cached-', myContext, '-', myAnnot)
  se = NULL
  if(dir.exists(dname)) {
    se = HDF5Array::loadHDF5SummarizedExperiment(dir=dname)
    flog.info('Read cached SE %s', dname)
  } else {
    flog.warn('Cached data for %s %s not found at %s, SKIP', myContext, myAnnot, dname)
    next
  }
  gpc_se <<- se
  ai = as.data.table(rowRanges(se))
  ai[['idx']] = seq(nrow(ai))
  gpcAnnoInfo <<- ai

  # Ensure column order
  if(!all(colnames(rnaNormData) == colnames(gpc_se))) {
    flog.warn('Fix RNA column order from %s', myAnnot)
    stopifnot(ncol(rnaNormData) <= ncol(gpc_se) && all(colnames(rnaNormData) %in% colnames(gpc_se)))
    nn = match(colnames(gpc_se), colnames(rnaNormData))
    nn = nn[!is.na(nn)]
    rnaNormData = rnaNormData[, nn]
    rnaCellInfo = rnaCellInfo[nn, ]
    stopifnot(all(colnames(rnaNormData) == colnames(gpc_se)))
  }
  
  if(featureData[[myAnnot]][['match']] == 'byPos') {
    
    opairs = getFeaturePairs(gpcAnnoInfo, rnaGeneInfo, stranded=TRUE, distance=featureData[[myAnnot]][['distance']])
    flog.info('Matching annotation %s by unstranded distance: %d pairs', myAnnot, nrow(opairs))
    
  } else if (featureData[[myAnnot]][['match']] == 'byID') {
    
    opairs = gpcAnnoInfo[rnaGeneInfo, on='ID==ensembl', nomatch=NULL][, .(xid=idx, yid=i.idx)]
    opairs = unique(opairs)
    flog.info('Matching annotation %s by ID: there are %d pairs for %d genes', myAnnot, nrow(opairs), nrow(rnaGeneInfo))
    
  }
  # Add a filter by minimum number of total reads per window per cell, derived from the total count matrix
  tMtx = assay(gpc_se, 'total')
  mMtx = assay(gpc_se, myMethAssay)
  nNA.i = sum(!is.na(mMtx))
  mMtx[tMtx < myWinTotalRC] = NA
  nNA.f = sum(!is.na(mMtx))
  flog.info('Minimum window read count filter of > %d total reads reduces the non-NA elements from %d to %d', myWinTotalRC, nNA.i, nNA.f)
  
  nCorrel = .pairedCorrel(aMtx=      mMtx,
                          bMtx=      rnaNormData,
                          pairing=   opairs,
                          groups=    unlist(rnaCellInfo[['Treatment']]),
                          minObsProp=mySampProp,
                          doBiweight=FALSE,
                          BPP=       BPP)
  aCorrel = cbind(gpcAnnoInfo[nCorrel[['aIdx']], .(anno.win = paste0(seqnames, ':', start, '-', end, '.(', strand, ')'),
                                        ID.x = ID)],
                  rnaGeneInfo[nCorrel[['bIdx']], .(ID.y = ensembl)],
                  nCorrel,
                  gpcAnnoInfo[nCorrel[['aIdx']], .(seqnames, start.x = start, end.x = end, strand.x = strand)],
                  rnaGeneInfo[nCorrel[['bIdx']], .(start.y = start, end.y = end, strand.y = strand, hgnc, biotype)])
  
  aKey = paste0(myContext, '-RNA:', myAnnot)
  corResults[[aKey]] = list(aVals=myMethAssay, bVals='RNANormNA', cor=aCorrel, methMinRC=myWinTotalRC,
                            match=featureData[[myAnnot]][['match']], distance=featureData[[myAnnot]][['distance']])
  # Save at each annotation step, in case something interrupts the process in between
  saveRDS(corResults, saveFile)
  
  flog.info('End process of %s: %d correlations on all groups', myAnnot, nrow(aCorrel))
  
}


stop(4)
