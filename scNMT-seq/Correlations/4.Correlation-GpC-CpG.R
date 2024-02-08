#
# Perform the methylation correlation analysis between methylation contexts
#
# We use the results previously computed and agreed:
# CpG data on 222 cells, rates not batch corrected.
# GpC data on 222 cells, normalised-batch corrected rates
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

source('methUtils.R')
source('./methCorUtils.R')

logFile = strftime(Sys.Date(), "correlation-%Y%m%d.log")
flog.appender(flog.indent_wrapper(appender=appender.tee(logFile), indent="-       "))
flog.info('==========================\n                  ** Start Correlation Analysis **\n                  ==========================\n')

mvpar = methRegionsPar(logDir=workDir, numWorkers=18, jobname='Corr')

mySampProp = 0.05
matchPolicy = 'byID'
myWinGTotalRC = 20
myWinCTotalRC = 5

mxplDir = paste0(workDir3, 'methXplor/')

BPP = getBPP(params=mvpar)

#############
# CpG - GpC correlation
flog.info('%%\n%%       Process GpC-CpG Context\n%%')

myContext = NULL
saveFile = paste0('correlations-GpC-CpG-3d.rds')

myAnnot = 'CGI_SeqMonk'
myAnnot = Sys.getenv('CORR_ANNOT', unset=NA)
if (!is.na(myAnnot)) {
  featureData = featureData[myAnnot]
}
corResults = list()
for (myAnnot in names(featureData)) {
  
  flog.info('\nAnnotation: %s\n', myAnnot)
  
  dname = paste0(mxplDir, 'cached-CpG-', myAnnot)
  seC = NULL
  if(dir.exists(dname)) {
    seC = HDF5Array::loadHDF5SummarizedExperiment(dir=dname)
    flog.info('Read cached SE %s', dname)
  } else {
    flog.warn('Cached data for CpG %s not found at %s, SKIP', myAnnot, dname)
    next
  }
  
  dname = paste0(mxplDir, 'cached-GpC-', myAnnot)
  seG = NULL
  if(dir.exists(dname)) {
    seG = HDF5Array::loadHDF5SummarizedExperiment(dir=dname)
    flog.info('Read cached SE %s', dname)
  } else {
    flog.warn('Cached data for GpC %s not found at %s, SKIP', myAnnot, dname)
    next
  }
  
  cpg_se <<- seC
  gpc_se <<- seG
  aiC = as.data.table(rowRanges(cpg_se))
  aiC[['idx']] = seq(nrow(aiC))
  cpgAnnoInfo <<- aiC
  aiG = as.data.table(rowRanges(gpc_se))
  aiG[['idx']] = seq(nrow(aiG))
  gpcAnnoInfo <<- aiG
  
  # Ensure column order
  if(!all(colnames(cpg_se) == colnames(gpc_se))) {
    flog.warn('Fix CpG column order from %s', myAnnot)
    stopifnot(ncol(cpg_se) <= ncol(gpc_se) && all(colnames(cpg_se) %in% colnames(gpc_se)))
    nn = match(colnames(gpc_se), colnames(cpg_se))
    nn = nn[!is.na(nn)]
    cpg_se = cpg_se[, nn]
    stopifnot(all(colnames(cpg_se) == colnames(gpc_se)))
  }
  
  # Ensure we have groups
  if (!is.null(colData(gpc_se)[['Treatment']]) && 
      !is.null(colData(cpg_se)[['Treatment']]) && 
      !all(colData(gpc_se)[['Treatment']] == colData(cpg_se)[['Treatment']])) {
    flog.error('Fatal error in sample labels')
    stop(5)
  } else if (is.null(colData(gpc_se)[['Treatment']]) && is.null(colData(cpg_se)[['Treatment']])) {
    flog.error('Missing group info !')
    stop(6)
  } else if (is.null(colData(gpc_se)[['Treatment']]) && !is.null(colData(cpg_se)[['Treatment']])) {
    colData(gpc_se)[['Treatment']] = colData(cpg_se)[['Treatment']]
  } else if (!is.null(colData(gpc_se)[['Treatment']]) && is.null(colData(cpg_se)[['Treatment']])) {
    colData(cpg_se)[['Treatment']] = colData(gpc_se)[['Treatment']]
  }
  stopifnot(!is.null(colData(gpc_se)[['Treatment']]) && all(colData(gpc_se)[['Treatment']] %in% c('AZA', 'DAC', 'Unt')))
  
  # CpG - GpC matches by ID
  if(matchPolicy == 'byDistance') {
    
    opairs = getFeaturePairs(gpcAnnoInfo, cpgAnnoInfo, stranded=FALSE, distance=featureData[[myAnnot]][['distance']])
    flog.info('Matching annotation %s by unstranded distance: %d pairs', myAnnot, nrow(opairs))
    
  } else if (matchPolicy == 'byID') {
    
    opairs = gpcAnnoInfo[cpgAnnoInfo, on=c('seqnames','start','end','strand'), nomatch=NULL][, .(xid=idx, yid=i.idx)]
    opairs = unique(opairs)
    flog.info('Matching annotation %s by strict equality: %d pairs', myAnnot, nrow(opairs))
    
  } else {
    flog.error("Unknown match policy: %s", matchPolicy)
    stop(3)
  }
  
  # Add a filter by minimum number of total reads per window per cell, derived from the total count matrix
  tMtx = assay(gpc_se, 'total')
  gMtx = assay(gpc_se, 'norm.rate')
  gnNA.i = sum(!is.na(gMtx))
  gMtx[tMtx < myWinGTotalRC] = NA
  gnNA.f = sum(!is.na(gMtx))
  flog.info('Minimum window read count filter of > %d total reads reduces the GpC non-NA elements from %d to %d', myWinGTotalRC, gnNA.i, gnNA.f)
  
  tMtx = assay(cpg_se, 'total')
  cMtx = assay(cpg_se, 'rate')
  cnNA.i = sum(!is.na(cMtx))
  cMtx[tMtx < myWinCTotalRC] = NA
  cnNA.f = sum(!is.na(cMtx))
  flog.info('Minimum window read count filter of > %d total reads reduces the CpG non-NA elements from %d to %d', myWinCTotalRC, cnNA.i, cnNA.f)
  
  nCorrel = .pairedCorrel(aMtx=      gMtx,
                          bMtx=      cMtx,
                          pairing=   opairs,
                          groups=    unlist(colData(gpc_se)[['Treatment']]),
                          minObsProp=mySampProp,
                          doBiweight=FALSE,
                          BPP=       BPP)
  aCorrel = cbind(gpcAnnoInfo[nCorrel[['aIdx']], .(gpc.win = paste0(seqnames, ':', start, '-', end, '.(', strand, ')'),
                                        ID.x = ID)],
                  cpgAnnoInfo[nCorrel[['bIdx']], .(cpg.win = paste0(seqnames, ':', start, '-', end, '.(', strand, ')'),
                                        ID.y = ID)],
                  nCorrel,
                  gpcAnnoInfo[nCorrel[['aIdx']], .(seqnames, start.x = start, end.x = end, strand.x = strand)],
                  cpgAnnoInfo[nCorrel[['bIdx']], .(start.y = start, end.y = end, strand.y = strand)])
  
  aKey = paste0('GpC-CpG:', myAnnot)
  corResults[[aKey]] = list(aVals='norm.rate', bVals='rate', cor=aCorrel, methMinRC=c(gpc=myWinGTotalRC, cpg=myWinCTotalRC),
                            match=matchPolicy, distance=featureData[[myAnnot]][['distance']])
  # Save at each annotation step, in case something interrupts the process in between
  saveRDS(corResults, saveFile)
  
  flog.info('End process of %s: %d correlations on all groups', myAnnot, nrow(aCorrel))
  
}



stop(4)
