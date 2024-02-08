#
# Prepare data for use in the methylation explorer.
#
# This script loads individual cell coverage data and normalises the methylation counts to obtain normalised methylation rates.
# The rates / counts are then filtered by annotation regions.
#
# Distilled from analysis on CpG data. 
#
# CR 20210121, 20220825, 20221031, 20230308
#
# Ok. So lay out the normalisation strategy.
# After much reading and deliberation it seems apparent that previous deconvolution normalisation strategy is not adequate, or at least insufficient.
# In this script I will attempt the following steps:
# 1. Filtering cells by QC PASS status as provided. Let's not pollute our field.
# 2. I think the read counts are extremely low at the locus (cytosine) level to be adequately modeled / corrected, therefore I will
#    use the "agnostic" 3kb windows at 1kb steps as pooling factors.  Mainly because is one of the annotations of interest, and
#    secondly because by having a 2/3ds overlap between consecutive windows, they provide a strong smoothing across the genome.
# 3. Initial per-cell scaling factor. Compute per cell factors on total read count from coverage, such that all cells in all batches
#    have an initial comparable library size. This should be computed on CpG + GpC counts, rather than by context. Maybe I will.
# 4. Now on context. regression of batch + tratment on scaled methylated and unmethylated counts by separate.
# 5. Investigate MNN
#
# Implementation notes:
# 1. Use of the collected matrices 'methylVariance-GpC-Annotation-unnorm.csv.gz, to avoid reading all coverage files again and 
#    grouping by annotation window
#
source('/home/carlos/Projects/BIO19003-201707-Canepi-Lee/docs/20221031-CorrelationAnalysis-3Batches-v2/definitions.R')

setwd(workDir)

library(futile.logger)
library(data.table)
library(SummarizedExperiment)
library(Matrix)
library(HDF5Array)

source('./methCorUtils.R')

logFile = strftime(Sys.Date(), "fullnormalise-%Y%m%d.log")
flog.appender(flog.indent_wrapper(appender=appender.tee(logFile), indent="-       "))

.inferColData <- function(cnames, split='_', colord=c(2,3,1)) {
  dd = rbindlist(lapply(strsplit(cnames, split=split), as.list))[, ..colord]
  colnames(dd) = c('Treatment','Batch','Well')
  dd[Batch == '060320', Batch := 'b320']
  dd[Batch == '250620', Batch := 'b620']
  dd[Batch == '220331', Batch := 'b322']
  dd$Well = gsub('^([0-9]+)([A-Z])$', '\\2\\1', dd$Well)
  dd$Well = gsub('^([A-Z])([0-9]{1})$', '\\10\\2', dd$Well)
  return(dd)
}

# Not sure I need this, but just in case
x = 'All_loci'
if(is.null(featureData[[x]])) {
  dd = get_all_std_chr_regions()
  featureData[[x]] = list(dataFile=NULL, adjustStart=F, distance=1, match='byPos', data=dd, step.max.rows=1000)
  flog.info('Added and process %s\n   Regions: %d\n   Unique IDs: %d\n   Unique UIDs: %d\n   Step size: %d', 
            x, nrow(dd), uniqueN(dd$ID), uniqueN(dd$UID), featureData[[x]][['step.max.rows']])
}

x = 'MyWinM5S25'
if(is.null(featureData[[x]])) {
  dd = get_all_std_chr_regions(winLen=500000L, winStep=250000L)
  featureData[[x]] = list(dataFile=NULL, adjustStart=F, distance=1, match='byPos', data=dd, step.max.rows=1000)
  flog.info('Added and process %s\n   Regions: %d\n   Unique IDs: %d\n   Unique UIDs: %d\n   Step size: %d', 
            x, nrow(dd), uniqueN(dd$ID), uniqueN(dd$UID), featureData[[x]][['step.max.rows']])
}

qcFile = paste0(qcDir, 'Passing_QC_all_batches.txt')
qcCellList = fread(qcFile)
ss = apply(.inferColData(qcCellList[['Sample']]), 1, paste0, collapse='.')
set(qcCellList, j=c('Treatment','Batch','Well','Sample'), value=c(as.list(.inferColData(qcCellList[['Sample']])), list(ss)))

mvpar = methRegionsPar(logDir=workDir, numWorkers=5, jobname='Norm')
BPP = getBPP(params=mvpar)

#### Methylation data ######
myContext = 'CpG'
# for(myContext in c('CpG','GpC')) {
  
  flog.info("----\n\t\t   ---- START Preprocess for context %s\n\t\t   ----", myContext)
  
  for(x in names(featureData)) {
    featureData[[x]]$unnFile = paste0(workDir3, 'methylVariance-', x, '-', myContext, '-unnorm.csv.gz')         # Per region counts
    featureData[[x]]$norFile = paste0(workDir3, 'methXplor/cached-', myContext, '-', x, '-norm.rate.csv.gz')    # Normalised rates
    featureData[[x]]$coldFile = paste0(workDir3, 'methXplor/cached-', myContext, '-', x, '-coldata.csv.gz')     # Column data
    featureData[[x]]$rowdFile = paste0(workDir3, 'methXplor/cached-', myContext, '-', x, '-rowdata.csv.gz')     # Column data
    featureData[[x]]$hdf5Dir = paste0(workDir, 'se-meth-', myContext, '-', x)
  }
  
  myFeature = 'MyWinM5S25'
  FD = featureData[[myFeature]]  
  
  if(!dir.exists(FD[['hdf5Dir']])) {
    
    if (!file.exists(FD[['unnFile']])) {
      flog.info('Unnormalised file for %s not found, regenerating (%s)', myFeature, FD[['unnFile']])
      covFiles = list.files(path=covDir, pattern=paste0(".*.NOMe.", myContext, ".cov.gz"), full.names=T, recursive=F)
      qcFF = .inferColData(sapply(basename(covFiles), function(z) unlist(strsplit(z, '_cytosine'))[1], USE.NAMES=F))
      qcnn = qcFF[qcCellList, on=c('Treatment','Batch','Well'), which=T]
      # This is just temporary
      qcCovFiles = covFiles[qcnn]
      
      # covFiles = sample(qcCovFiles, 10)
      # qcFF = .inferColData(sapply(basename(covFiles), function(z) unlist(strsplit(z, '_cytosine'))[1], USE.NAMES=F))
      
      ccd = apply(qcFF, 1, paste0, collapse='.')

      methRegions = FD$data
      methVar = methVar_dataset(covFiles, methRegions, export_all=T, params=mvpar)
      flog.info('From coverage: %d regions', nrow(methVar))
      cns = colnames(methVar)
      cns = cns[(ncol(methRegions)):length(cns)]
      cnw = sapply(cns, function(t) unlist(strsplit(t,'.',fixed=T))[1], USE.NAMES=F)
      ccd = apply(.inferColData(cnames=cnw), 1, paste0, collapse='.')
      cww = sapply(seq(length(ccd)), function(k) gsub(paste0(cnw[k], '.'), paste0(ccd[k], '_'), cns[k], fixed=T))
      setnames(methVar, old=cns, new=cww)
      nn = which(apply(methVar[, .SD, .SDcols=grep('total.n', colnames(methVar))], 1, function(z) { sum(is.na(z)) >= 22 }))
      if(length(nn) > 0) {
        flog.info('Removing %d regions with more than 10%% of missing values', length(nn))
        flog.info('', methVar[nn, UID], capture=T)
        methVar = methVar[-nn, ]
      }
      fwrite(methVar, file=FD[['unnFile']])
      rm(methVar)
    }
    
    FD$exps = SumExpFromRateFile(rateFile=FD[['unnFile']], joinchar='_')
    # Just those with QC PASS
    FD$exps = FD$exps[, qcCellList[['Sample']]]
    cd = .inferColData(colnames(FD$exps), split='[.]', colord=c(1,2,3))
    dd = as.matrix(assay(FD$exps, i='total'))
    nn = which(apply(dd, 1, function(a) all(is.na(a))))
    if(length(nn) > 0) {
      dd = dd[-nn, ]
      FD$exps = FD$exps[-nn, ]
    }
    set(cd, j=c('nRgns','totLociCover'), value=list(colSums(!is.na(dd)), colSums(dd, na.rm=T)))
    dd = as.matrix(assay(FD$exps, i='loci'))
    set(cd, j=c('nLoci'), value=list(colSums(dd, na.rm=T)))
    colData(FD$exps) = cbind(colData(FD$exps), as(cd, 'DataFrame'))
    
    # Standardise storage modes
    for (nn in assayNames(FD$exps)) {
      aa = assay(FD$exps, nn)
      if(is.data.table(aa)) {
        aa = Matrix::Matrix(as.matrix(aa))
        assay(FD$exps, nn, withDimnames=F) = aa
      }
    }
    # Save object in quicker access format
    saveHDF5SummarizedExperiment(FD$exps, dir=FD[['hdf5Dir']], verbose=T)
    rm(FD)
    FD = featureData[[myFeature]]  
  }
  
  FD$exps = loadHDF5SummarizedExperiment(dir=FD[['hdf5Dir']])
  
  if(!all(c('met.norm','unmet.norm') %in% assayNames(FD$exps))) {
    cd = as.data.table(colData(FD$exps))
    flog.info('Read region methylation info for "%s" from: %s', myFeature, basename(FD[['hdf5Dir']]))
    flog.info('%d rows by %d columns', nrow(FD$exps), ncol(FD$exps))
    treatClusters = unlist(cd[ , .(ifelse(Treatment %in% c('Neg','Unt'), 'No', Treatment))])
    
    doRep <- function(ilist, pdd) {
      ir = ic = NULL
      for(jc in seq(ncol(pdd))) {
        i = ilist - (jc-1) * nrow(pdd)
        i = i[between(i, 1, nrow(pdd))]
        if(length(i) == 0) next
        ir = c(ir, i)
        ic = c(ic, rep(jc, length(i)))
        # ddd[i, jc] = pdd[i, jc]
      }
      return(list(r=ir, c=ic))
    }
    
    # Not well behaved on NAs. We have a small number of NAs and will impute by kNN.  Imputed value is written back to assay
    dd = assay(FD$exps, i='total')
    if(sum(is.na(dd)) > 0) {
      flog.info('Total counts: %d NAs, %d non-NAs', sum(is.na(dd)), sum(!is.na(dd)))
      idd = impute::impute.knn(as.matrix(dd), rng.seed=124986709)
      storage.mode(idd$data) = type(dd)
      ix = doRep(which(is.na(dd)), idd$data)
      # Complicated way to write back to not loose the connection with original HDF5 array
      dd = assay(FD$exps, 'total')
      dd[ix$r, ix$c] = idd$data[ix$r, ix$c]
      assay(FD$exps, 'total') = dd
    } else {
      flog.info('No NAs in "total"')
    }
    qq = scran::quickCluster(dd, min.size=12, d=NA, use.ranks=T, block=cd$Batch, BPPARAM=BPP)
    flog.info('QuickClusters with batch blocking returns %d clusters', length(levels(qq)))
    flog.info('Clusters vs. Batch, Treatment', table(Clus=qq, paste(cd$Batch, cd$Treatment)), capture=T)
    FD$exps$qCluster = qq
    FD$exps$sf.lib = scuttle::librarySizeFactors(dd, BPPARAM=BPP)
    FD$exps$sf.geom = scuttle::geometricSizeFactors(dd, BPPARAM=BPP)
    
    # We will correct for batch later on, so our clusters are at least batch, or subgroups in batch.
    # Reduce the sizes from default values due to the distribution of treatments, 
    # eg check table(cd[, c('Treatment','Batch')]), table(cd[, c('Batch','qCluster')])
    FD$exps$sf.pool0   = scuttle::pooledSizeFactors(dd, sizes=seq(21,101,5), BPPARAM=BPP)
    FD$exps$sf.poolB   = scuttle::pooledSizeFactors(dd, sizes=seq(11,91,5), cluster=cd$Batch,    BPPARAM=BPP)
    FD$exps$sf.poolQ   = scuttle::pooledSizeFactors(dd, sizes=seq(11,91,5), cluster=cd$qCluster, BPPARAM=BPP)
    FD$exps$sf.poolB2  = scuttle::pooledSizeFactors(dd, sizes=seq(11,91,5), scaling=FD$exps$sf.pool0, BPPARAM=BPP)
    flog.info('%s: %s Size factors on total counts', myFeature, myContext)
    
    dd = assay(FD$exps, i='meth')
    if(sum(is.na(dd)) > 0) {
      flog.info('Meth counts: %d NAs, %d non-NAs', sum(is.na(dd)), sum(!is.na(dd)))
      idd = impute::impute.knn(as.matrix(dd), rng.seed=908123642)
      storage.mode(idd$data) = type(dd)
      ix = doRep(which(is.na(dd)), idd$data)
      # Complicated way to write back to not loose the connection with original HDF5 array
      dd = assay(FD$exps, 'meth')
      dd[ix$r, ix$c] = idd$data[ix$r, ix$c]
      assay(FD$exps, 'meth') = dd
    } else {
      flog.info('No NAs in "meth"')
    }
    FD$exps$sfm.lib     = scuttle::librarySizeFactors(dd, BPPARAM=BPP)
    FD$exps$sfm.pool0   = scuttle::pooledSizeFactors(dd, sizes=seq(21,101,5), BPPARAM=BPP)
    FD$exps$sfm.poolS   = scuttle::pooledSizeFactors(dd, sizes=seq(21,101,5), scaling=FD$exps$sfm.lib, BPPARAM=BPP)
    FD$exps$sfm.poolB   = scuttle::pooledSizeFactors(dd, sizes=seq(11,91,5), cluster=cd$Batch, BPPARAM=BPP)
    FD$exps$sfm.poolBS  = scuttle::pooledSizeFactors(dd, sizes=seq(11,91,5), cluster=cd$Batch, scaling=FD$exps$sf.lib, BPPARAM=BPP)
    FD$exps$sfm.poolBS2 = scuttle::pooledSizeFactors(dd, sizes=seq(11,91,5), cluster=cd$Batch, scaling=FD$exps$sfm.poolBS, BPPARAM=BPP)
    flog.info('%s: %s Size factors on methylated counts', myFeature, myContext)
    
    dd = as.matrix(assay(FD$exps, i='unmeth'))
    if(sum(is.na(dd)) > 0) {
      flog.info('UnMeth counts: %d NAs, %d non-NAs', sum(is.na(dd)), sum(!is.na(dd)))
      idd = impute::impute.knn(as.matrix(dd), rng.seed=23984)
      storage.mode(idd$data) = type(dd)
      ix = doRep(which(is.na(dd)), idd$data)
      # Complicated way to write back to not loose the connection with original HDF5 array
      dd = assay(FD$exps, 'unmeth')
      dd[ix$r, ix$c] = idd$data[ix$r, ix$c]
      assay(FD$exps, 'unmeth') = dd
    } else {
      flog.info('No NAs in "unmeth"')
    }
    FD$exps$sfu.lib     = scuttle::librarySizeFactors(dd, BPPARAM=BPP)
    FD$exps$sfu.pool0   = scuttle::pooledSizeFactors(dd, sizes=seq(21,101,5), BPPARAM=BPP)
    FD$exps$sfu.poolS   = scuttle::pooledSizeFactors(dd, sizes=seq(21,101,5), scaling=FD$exps$sfu.lib, BPPARAM=BPP)
    FD$exps$sfu.poolB   = scuttle::pooledSizeFactors(dd, sizes=seq(11,91,5), cluster=cd$Batch, BPPARAM=BPP)
    FD$exps$sfu.poolBS  = scuttle::pooledSizeFactors(dd, sizes=seq(11,91,5), cluster=cd$Batch, scaling=FD$exps$sf.lib, BPPARAM=BPP)
    FD$exps$sfu.poolBS2 = scuttle::pooledSizeFactors(dd, sizes=seq(11,91,5), cluster=cd$Batch, scaling=FD$exps$sfm.poolBS, BPPARAM=BPP)
    flog.info('%s: %s Size factors on unmethylated counts', myFeature, myContext)
    
    rifiplot <- function(vx='sf.lib', vy='sf.geom', vf='Batch', cd=as.data.table(colData(FD$exps))) {
      require(ggplot2)
      gg = ggplot2::ggplot(data=cd) + ggplot2::geom_point(aes(x=.data[[vx]], y=.data[[vy]], fill=.data[[vf]], text=.data[['sample']]), stroke=.2, alpha=.6) + 
        ggplot2::geom_abline(slope=1, color='blue', size=.1) + theme_classic()
      plotly::ggplotly(gg)
    }
    
    dm = scuttle::normalizeCounts(assay(FD$exps, 'meth'),   transform='none', size_factors=FD$exps$sfm.poolBS, BPPARAM=BPP)
    du = scuttle::normalizeCounts(assay(FD$exps, 'unmeth'), transform='none', size_factors=FD$exps$sfu.poolBS, BPPARAM=BPP)
    
    assay(FD$exps, 'met.norm', withDimnames=F)   = dm
    assay(FD$exps, 'unmet.norm', withDimnames=F) = du
    
    quickResaveHDF5SummarizedExperiment(FD$exps, verbose=T)
  }
  
  if(!all(c('sfm.libN','sfm.bc1','sfu.libN','sfu.bc1') %in% colnames(colData(FD$exps)))) {
    # This seems to get stuck
    if(FALSE) {
      dm = assay(FD$exps, 'met.norm')
      dm[is.na(dm)] = 0
      snnmg = scran::buildSNNGraph(dm, 
                                   BNPARAM=BiocNeighbors::KmknnParam(distance='Manhattan'),
                                   BSPARAM=BiocSingular::RandomParam(), 
                                   BPPARAM=BPP)
      clus.m.norm = igraph::cluster_walktrap(snnmg)
      mixt.norm = table(clus.m.norm=clus.m.norm$membership, batch=colData(FD$exps)$Batch)
    }
    
    # Try rescaling batch on normalised counts
    dm = scuttle::normalizeCounts(assay(FD$exps, 'meth'), transform='log', size_factors=FD$exps$sfm.poolBS, BPPARAM=BPP)
    du = scuttle::normalizeCounts(assay(FD$exps, 'unmeth'), transform='log', size_factors=FD$exps$sfu.poolBS, BPPARAM=BPP)
    dmbc = batchelor::rescaleBatches(dm, batch=FD$exps$Batch)
    dubc = batchelor::rescaleBatches(du, batch=FD$exps$Batch)
    
    assay(FD$exps, 'meth.s.bc') = assay(dmbc, 'corrected')
    assay(FD$exps, 'unmeth.s.bc') = assay(dubc, 'corrected')
    
    dd = assay(FD$exps, 'met.norm')
    sfbc = colSums(dd, na.rm=T)
    sfbc = sfbc / mean(sfbc)
    FD$exps$sfm.libN = sfbc
    
    dd = as.data.table(2^assay(dmbc, 'corrected'))
    sfbc = colSums(dd, na.rm=T)
    sfbc = sfbc / mean(sfbc)
    FD$exps$sfm.bc1 = sfbc
    
    dd = assay(FD$exps, 'unmet.norm')
    sfbc = colSums(dd, na.rm=T)
    sfbc = sfbc / mean(sfbc)
    FD$exps$sfu.libN = sfbc
    
    dd = as.data.table(2^assay(dubc, 'corrected'))
    sfbc = colSums(dd, na.rm=T)
    sfbc = sfbc / mean(sfbc)
    FD$exps$sfu.bc1 = sfbc
    
    # rifiplots
    
    quickResaveHDF5SummarizedExperiment(FD$exps, verbose=T)
  }
  
  if(!all(c('sfm.rgbc','sfu.rgbc') %in% colnames(colData(FD$exps)))) {

    dm = scuttle::normalizeCounts(assay(FD$exps, 'meth'), transform='log', size_factors=FD$exps$sfm.poolBS, BPPARAM=BPP)
    du = scuttle::normalizeCounts(assay(FD$exps, 'unmeth'), transform='log', size_factors=FD$exps$sfu.poolBS, BPPARAM=BPP)

    aa = colData(FD$exps)
    aa$Batch = factor(aa$Batch, levels=c('b322', 'b320', 'b620'))
    aa$Treatment = factor(aa$Treatment, levels=c('Unt', 'AZA', 'DAC'))
    aamm = model.matrix(~ Batch + Treatment, data=aa)
    
    drbc = batchelor::regressBatches(dm, batch=aa$Batch, design=aamm, keep=grep('Treatment|Intercept', colnames(aamm)), 
                                     d=NA, correct.all=T, deferred=F, BPPARAM=BPP)
    dmr = assay(drbc, 'corrected')
    dd = as.matrix(dm)
    dd[] = 0.0
    # Takes forever
    for(j in seq(nrow(dd))) {
      dd[j,] = dmr[j,]
      if(j %% 2000 == 1) cat(j, '..')
    }
    cat('\n')
    sfbc = colSums(2^dd)
    sfbc = sfbc / mean(sfbc)
    FD$exps$sfm.rgbc = sfbc
    
    drbc = batchelor::regressBatches(du, batch=aa$Batch, design=aamm, keep=grep('Treatment|Intercept', colnames(aamm)), 
                                     d=NA, correct.all=T, BPPARAM=BPP)
    dur = assay(drbc, 'corrected')
    dd = as.matrix(du)
    dd[] = 0.0
    for(j in seq(nrow(dd))) {
      dd[j,] = dur[j,]
      if(j %% 2000 == 1) cat(j, '..')
    }
    cat('\n')
    sfbc = colSums(2^dd)
    sfbc = sfbc / mean(sfbc)
    FD$exps$sfu.rgbc = sfbc
    
    quickResaveHDF5SummarizedExperiment(FD$exps, verbose=T)
    
    # takes forever
    if(FALSE) {
      # Normalised and batch corrected data
      dd = assay(drbc, 'corrected')
      set(dd, j='rn', value=list(rownames(FD$exps)))
      dz = melt.data.table(dd, id.vars='rn', na.rm=T)
      set(dz, j='Batch', value='b322')
      set(dz, i=grep('b320', dz$variable), j='Batch', value='b320')
      set(dz, i=grep('b620', dz$variable), j='Batch', value='b620')
      set(dz, j='Treat', value='DAC')
      set(dz, i=grep('AZA', dz$variable), j='Treat', value='AZA')
      set(dz, i=grep('Unt', dz$variable), j='Treat', value='Unt')
      set(dz, j='variable', value=NULL)
      gg = ggplot(data=dz) + 
        geom_boxplot(aes(x=Batch, y=value, fill=Batch), outlier.alpha=0.2, outlier.color='#6c757d') + 
        theme_classic() + 
        facet_wrap(~Treat, nrow=1)
      gg + ylim(1e-4, 15) + scale_y_log10()
    }
  }
  
  #Mutual Nearest neighbors batch correction
  # No idea how to use this as normalised data

  # Plot (takes an insane amount of time and memory)
  if(FALSE) {
    # Normalised, not batch corrected data, on same scale (log2)
    dd = as.data.table(du)
    set(dd, j='rn', value=list(rownames(FD$exps)))
    dz = melt.data.table(dd, id.vars='rn', na.rm=T)
    set(dz, j='Batch', value='b322')
    set(dz, i=grep('b320', dz$variable), j='Batch', value='b320')
    set(dz, i=grep('b620', dz$variable), j='Batch', value='b620')
    set(dz, j='Treat', value='DAC')
    set(dz, i=grep('AZA', dz$variable), j='Treat', value='AZA')
    set(dz, i=grep('Unt', dz$variable), j='Treat', value='Unt')
    set(dz, j='variable', value=NULL)
    gg = ggplot(data=dz) + 
      geom_boxplot(aes(x=Batch, y=value, fill=Batch), outlier.alpha=0.2, outlier.color='#6c757d') + 
      theme_classic() + 
      facet_wrap(~Treat, nrow=1)
    gg + ylim(1e-4, 15) + scale_y_log10()
  }
  if(FALSE) {
    # Normalised and batch corrected data
    dd = as.data.table(assay(dubc, 'corrected'))
    set(dd, j='rn', value=list(rownames(FD$exps)))
    dz = melt.data.table(dd, id.vars='rn', na.rm=T)
    set(dz, j='Batch', value='b322')
    set(dz, i=grep('b320', dz$variable), j='Batch', value='b320')
    set(dz, i=grep('b620', dz$variable), j='Batch', value='b620')
    set(dz, j='Treat', value='DAC')
    set(dz, i=grep('AZA', dz$variable), j='Treat', value='AZA')
    set(dz, i=grep('Unt', dz$variable), j='Treat', value='Unt')
    set(dz, j='variable', value=NULL)
    gg = ggplot(data=dz) + 
      geom_boxplot(aes(x=Batch, y=value, fill=Batch), outlier.alpha=0.2, outlier.color='#6c757d') + 
      theme_classic() + 
      facet_wrap(~Treat, nrow=1)
    gg + ylim(1e-4, 15) + scale_y_log10()
  }

  # Compute now rates and variance on scaled and scaled + bc counts
  # Use the normalised values
  dm = assay(FD$exps, 'met.norm')
  du = assay(FD$exps, 'unmet.norm')
  dr = (dm + 1) / (dm + du + 2)
  dv = dr * (1 - dr) / (dm + du)
  assay(FD$exps, 'rate.N') = dr
  assay(FD$exps, 'rvar.N') = dv
  
  dm = assay(dmbc, 'corrected')
  du = assay(dubc, 'corrected')
  dm = 2^dm - 1
  du = 2^du - 1
  drc = (dm + 1) / (dm + du + 2)
  dvc = dr * (1 - dr) / (dm + du)
  assay(FD$exps, 'rate.NBC') = drc
  assay(FD$exps, 'rvar.NBC') = dvc
  
  quickResaveHDF5SummarizedExperiment(FD$exps, verbose=T)
 
  # Thsi completes the estimation of the normalisation and batch correction factors for GpC.
  # Now we nee to apply to annotations and get the corrected rates and variances, generate the annotated datasets
  
#}

# At the end of the loop, the following quantities and size factors are calculated per cell:
# "sample", "idx", "Treatment", "Batch", "Well"     Meaning as usual
# "nRgns"         Number of regions in `MyWinM5S25` with data
# "totLociCover"  Accumulated read cover over loci in region
# "nLoci"         Number of loci with data in regions 
# "qCluster"      Cluster assignment from scran::QuickCluster()
# "sf.lib"        Library size factor on total region counts (m+u) in cell
# "sf.geom"       Geometric size factors on total region counts scuttle::geometricSizeFactors(). Not used in any computation
# "sf.pool0"      Pooled size factors (deconvolution SFs) on total region counts, no clusters or initial scaling.
# "sf.poolB"      Pooled size factors on total region counts, clustered by Batch.
# "sf.poolQ"      Pooled size factors on total region counts, clustered by quickClusters
# "sf.poolB2"     Pooled size factors on total region counts, initial scaling by sf.pool0
# "sf{m,u}.lib"      Library size factor on meth / unmeth region counts
# "sf{m,u}.pool0"    Pooled size factor on meth / unmeth region counts, no clustering or scaling
# "sf{m,u}.poolS"    Pooled size factor on meth / unmeth region counts, scaling by sf{m.u}.lib
# "sf{m,u}.poolB"    Pooled size factor on meth / unmeth region counts, cluster by Batch
# "sf{m,u}.poolBS"   Pooled size factor on meth / unmeth region counts, cluster by Batch, scaling by sf.lib.  USED FOR NORMALISATION
# "sf{m,u}.poolBS2"  Pooled size factor on meth / unmeth region counts, cluster by Batch, scaling by sf{m,u}.poolBS
# "sf{m,u}.bc1"      Scale factor from batchelor::rescaleBatches() on meth / unmeth region counts normalised by sf{m,u}.poolBS
# "sf{m,u}.libN"     Library size factor on meth / unmeth region counts normalised by sf{m,u}.poolBS
# "sf{m,u}.rgbc"     Scale factor from batchelor::regressBatches() after correction for Batch, computed to compare to sf{m.u}.bc1
#
# For the batchelor::regressBatches() size factors, the regression model is ~ Batch + Treatment, with b322 and Unt as reference levels,
# and the residual correction keeps Treatment and Intercept terms (that is, removes the differential Batch effects wrt. the reference level)
# The result is parallel to the sf{m.u}.bc1 factors, with the regression factors being smaller for b322 than the rescaledBatch ones
# (The opposite in b320)
#
# The assays in MyWinM5S25 are:
# loci       Number of loci covered per window (row)
# meth       Methylated counts in window
# unmeth     Unmethylated counts in window
# total      M+U
# rate       Raw rate in the binomial approximation
# var        Raw variance
# met.norm   Normalised meth log2 counts with sfm.poolBS
# unmet.norm Normalised unmethylated log2 counts with sfu.poolBS
# meth.s.bc  Methylated scaled and batch corrected (rescaleBatches) pseudo-counts
# unmeth.s.bc Unmethylated scaled and batch corrected (rescaleBatches) pseudo-counts
# rate.N     Normalised rate estimate
# rvar.N     Normalised variance estimate
# rate.NBC   Normalised and batch corrected estimate
# rvar.NBC   Normalised and batch corrected variance estimate

