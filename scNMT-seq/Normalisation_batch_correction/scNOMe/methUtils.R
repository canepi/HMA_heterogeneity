# Utility functions for methylation data
#
# CR, 20210212
#


# Function to do a panel plotting from a list of SummarizedExperiments with .coltypes assays.
quick.adhoc.plot <- function(featureData, dclass='total') {
  aSmry = data.table()
  aOutl = data.table()
  for(k in seq_along(featureData)) {
    cc = colData(featureData[[k]]$exps)
    batch = gsub('b..(0[36])..', 'b\\1', cc$Batch)
    group = paste0(cc$Group, '.', batch)
    dd = as.data.table(assay(featureData[[k]]$exps, dclass))
    for(g in unique(group)) {
      nn = which(colnames(dd) %in% rownames(cc)[group == g])
      # groupdata = log10(unlist(dd[, .SD, .SDcols=nn], use.names=F))
      groupdata = unlist(dd[, .SD, .SDcols=nn], use.names=F)
      ff = fivenum(groupdata)
      # Compute box & whisker fences
      iqr = ff[4]-ff[2]+1
      drange = c(max(ff[1], ff[2]-1.5*iqr), min(ff[5], ff[4]+1.5*iqr))
      no = which(!between(groupdata, drange[1], drange[2]))
      ww = mean(groupdata, na.rm=T)
      # Adjust max and min
      if(length(no) > 0) {
        ff[c(1,5)] = range(groupdata[-no], na.rm=T)
        aOutl = rbind(aOutl, data.table(feature=names(featureData)[k], group=g, data=groupdata[no]))
      }
      aSmry = rbind(aSmry, data.table(feature=names(featureData)[k], group=g, 
                                      min=ff[1], q1=ff[2], med=ff[3], q3=ff[4], max=ff[5], mean=ww))
    }
  }
  
  aSmry[, ':='(treat=substr(aSmry$group, 1, 3), batch=substr(aSmry$group, 5, 7))]
  gg = ggplot(data=aSmry, aes(x=group)) + 
    geom_boxplot(aes(ymin=min, lower=q1, middle=med, upper=q3, ymax=max, color=treat),
                 stat='identity')
  
  if(nrow(aOutl) > 0) {
    setorder(aOutl, feature, group, -data)
    aOutl[, ':='(treat=substr(group, 1, 3), batch=substr(group, 5, 7), i=.I)]
    nn = NULL
    for(p in unique(aOutl$feature)) {
      for(q in unique(aOutl$group)) {
        kk = aOutl[feature == p & group == q, i]
        if(length(kk) > 500) {
          nn = c(nn, sample(kk, 500), kk[1])   # Keep the first, which is maximum
        }
      }
    }
    if(length(nn) > 0) {
      aOutl = aOutl[nn,]
    }
    
    gg = gg + geom_point(data=aOutl, aes(x=group, y=data, group=batch), alpha=0.1, size=0.5)
  }
  
  gg = gg + ggplot2::facet_wrap(~feature, ncol=4)
  if(!grepl('rate', dclass)) {
    gg = gg + scale_y_log10()
  }
  title = paste0(dclass, ' counts')
  ylab = 'Counts'
  if(dclass=='rate') {
    title = dclass
    ylab = 'Meth. Rate'
  } else if(dclass=='var') {
    title = 'rate variance'
    ylab = 'Meth. Rate Variance'
  }
  gg = gg +
    apatheme + theme(axis.text.x=element_text(angle=90)) + 
    ggtitle(title) + ylab(ylab)
  gg
}

# Stripped down version of sampleqc data
.getSQC <- function(dataDir) {
  qcspreadsheetFile = paste0(dataDir, 'scNMTseq_QC.xlsx')
  
  qcd1 = data.table::data.table(openxlsx::readWorkbook(xlsxFile=qcspreadsheetFile, sheet='RNA_Seq_QC'))
  colnames(qcd1) = paste0('RNA.', make.names(colnames(qcd1)))
  invisible(data.table::setalloccol(qcd1))
  sns1 = gsub('_March[_S0-9]*$', '_060320', qcd1$RNA.Sample)
  sns1 = gsub('_June[_S0-9]*$', '_250620', sns1)
  data.table::set(qcd1, j=c('Sample', 'RNA.X12', 'RNA.X13'), value=list(sns1, NULL, NULL))
  data.table::setorder(qcd1, Sample)
  
  qcd2 = data.table::data.table(openxlsx::readWorkbook(xlsxFile=qcspreadsheetFile, sheet='NOME_Seq_QC'))
  colnames(qcd2) = paste0('NOME.', make.names(colnames(qcd2)))
  invisible(data.table::setalloccol(qcd2))
  sns2 = gsub('_March[_S0-9]*$', '_060320', qcd2$NOME.Sample)
  sns2 = gsub('_June[_S0-9]*$', '_250620', sns2)
  sns2 = gsub('Untreated', 'Unt', sns2)
  sns2 = gsub('Negative', 'Neg', sns2)
  data.table::set(qcd2, j=c('Sample', 'NOME.X13'), value=list(sns2, NULL))
  data.table::setorder(qcd2, Sample)
  sampleQC = qcd1[qcd2, on='Sample']
  sampleQC[, Pass := F]
  sampleQC[, Pass := (RNA.PASS.FAIL=='PASS' & NOME.QC.status=='PASS')]
  sampleQC = sampleQC[Pass==T, ]
  sampleQC[, ':='(Group = as.character(lapply(Sample, function(z) unlist(strsplit(z,'_'))[2])),
                  Batch = lapply(sampleQC$Sample, function(x) gsub('.*_[0-9]{2}([0-9]{2})[0-9]{2}$', 'b\\1', x)))]
  rownames(sampleQC) = sampleQC$Sample
  return(sampleQC)
}


#' computePairedCorrelation
#' 
#' @description Computes the pairwise correlation test on two matrices A and B indexed by am index pairing mapping.
#' 
#' @param matrixA,matrixB  Matrices to be correlated pairwise by paired indexed rows
#' @param pairing          A two-column data frame or matrix containing the row index pairs for matrix A and B
#' @param weightA,weightB  Optional weight matrices for weighted correlation analysis. If given must be same dimension than 
#'                         corresponding A or B matrix
#' @param use   The \code{use} parameter of WGCNA::corAndPvalue()
#' @param method  Correlation method to use.
#' @param BPP     An optional BiocParallel::BiocParallelParam() object for parallel computation.  If NULL a MulticoreParam with 15 workers is used
#' 
computePairedCorrelation <- function(matrixA, matrixB, pairing, weightA=NULL, weightB=NULL, 
                                     use='pairwise.complete.obs', method='pearson', BPP=NULL) {
  pairing = as.data.frame(pairing)
  if(ncol(matrixA) != ncol(matrixB) 
     || (!is.null(weightA) && ncol(weightA) != ncol(matrixA))
     || (!is.null(weightB) && ncol(weightB) != ncol(matrixB))) {
    stop('Not all data or weight matrices have same number of columns')
  } else if((!is.null(weightA) && nrow(weightA) != nrow(matrixA))
            || (!is.null(weightB) && nrow(weightB) != nrow(matrixB))) {
    stop('Weights and data have different number of rows')
  } else if(ncol(pairing) < 2 || !is.numeric(pairing[,1]) || !is.numeric(pairing[,2]) 
            || nrow(matrixA) < max(pairing[,1]) || nrow(matrixB) < max(pairing[,2])) {
    stop('Data rows and pairing index disagree')
  }
  
  if(is.null(BPP)) {
    BPP = BiocParallel::MulticoreParam(workers=15, progressbar=F)
  }
  
  # Correlation cases:
  if(!is.null(weightA) && !is.null(weightB)) {
    # Both weights non-null
    myCorrel = data.table::rbindlist(BiocParallel::bplapply(seq(nrow(pairing)), function(i, p, A, B, wA, wB, method, use) {
      aa = A[p[i,1], ] * wA[p[i,1], ]
      bb = B[p[i,2], ] * wB[p[i,2], ]
      WGCNA::corAndPvalue(aa, bb, method=method, use=use)
    }, p=pairing, A=matrixA, B=matrixB, wA=weightA, wB=weightB, method=method, use=use, BPPARAM=BPP))
  } else if(!is.null(weightA)) {
    # weight A non-null
    myCorrel = data.table::rbindlist(BiocParallel::bplapply(seq(nrow(pairing)), function(i, p, A, B, ww, method, use) {
      aa = A[p[i,1], ] * ww[p[i,1], ]
      WGCNA::corAndPvalue(aa, B[p[i,2], ], method=method, use=use)
    }, p=pairing, A=matrixA, B=matrixB, ww=weightA, method=method, use=use, BPPARAM=BPP))
  } else if(!is.null(weightB)) {
    # weight B non-null
    myCorrel = data.table::rbindlist(BiocParallel::bplapply(seq(nrow(pairing)), function(i, p, A, B, ww, method, use) {
      bb = B[p[i,2], ] * ww[p[i,2], ]
      WGCNA::corAndPvalue(A[p[i,1], ], bb, method=method, use=use)
    }, p=pairing, A=matrixA, B=matrixB, ww=weightB, method=method, use=use, BPPARAM=BPP))
  } else {
    # plain case
    myCorrel = data.table::rbindlist(BiocParallel::bplapply(seq(nrow(pairing)), function(i, p, A, B, method, use) {
      WGCNA::corAndPvalue(A[p[i,1], ], B[p[i,2], ], method=method, use=use)
    }, p=pairing, A=matrixA, B=matrixB, method=method, use=use, BPPARAM=BPP))
  }
  data.table::set(myCorrel, j=c('p.adj'), value=p.adjust(myCorrel$p, method='BH'))
  
  return(myCorrel)
}

#' checkPairingCompleteObs
#' 
#' Returns the indices of the pairing table that produce at least minobs non NA values from the combined matrices A and B
#' 
checkPairingCompleteObs <- function(matrixA, matrixB, pairing, minobs=3) {
  pairing = as.data.frame(pairing)
  if(ncol(matrixA) != ncol(matrixB)) {
    stop('Data matrices have different number of columns')
  } else if(ncol(pairing) < 2 || !is.numeric(pairing[,1]) || !is.numeric(pairing[,2]) 
            || nrow(matrixA) < max(pairing[,1]) || nrow(matrixB) < max(pairing[,2])) {
    stop('Data rows and pairing index disagree')
  }
  
  aMtx = Matrix::Matrix(is.na(matrixA)[pairing[, 1], ], nrow=nrow(pairing), sparse=T)
  bMtx = Matrix::Matrix(is.na(matrixB)[pairing[, 2], ], nrow=nrow(pairing), sparse=T)
  cMtx = !(aMtx | bMtx)
  nn = rowSums(cMtx)
  nn = which(nn >= minobs)
  return(nn)
}

# Auxiliary function, most data from the enclosing context
# Columns in input data and groups MUST be aligned, e.g. in same (column) order
# If groups is not NULL, the correlation will ALSO be performed in column subsets of the matrices.  
# In such case, results will be reported for all rows that fullfill the minObsProp simultaneously for all groups
.pairedCorrel <- function(aMtx, bMtx, pairing, groups=NULL, minObsProp=0.1, doBiweight=T, BPP=BiocParallel::bpparam()) {
  
  if(ncol(aMtx) != ncol(bMtx) || max(unlist(pairing[, 1])) > nrow(aMtx) || max(unlist(pairing[, 2])) > nrow(bMtx)) {
    stop('Wrong matrix dimensions or bad pairing')
  }
  cases = list(list(ccn=seq(ncol(aMtx)), label='', minobs=max(3, round(ncol(aMtx) * minObsProp))))
  if(!is.null(groups)) {
    groups = as.factor(groups)
    groups = droplevels.factor(groups)
    if(length(groups) != ncol(aMtx)) {
      stop('Groups length different from columns number')
    }
    cases = c(cases,
              lapply(levels(groups), function(x, groups) { 
                list(ccn=which(groups==x), label=as.character(x), minobs=max(3, round(sum(groups==x) * minObsProp))) },
                groups=groups)
    )
  }

  myCorrel = NULL
  BiocParallel::bpprogressbar(BPP) = T
  
  for(xl in cases) {  
    nnx = checkPairingCompleteObs(aMtx[, xl$ccn], bMtx[, xl$ccn], pairing=pairing, minobs=xl$minobs)
    flog.info('There are %d pairs with less than %d observations in group %s', nrow(pairing) - length(nnx), xl$minobs,
              ifelse(xl$label=='', '[ALL]', xl$label))
    flog.info('Computing correlation on %d pairs in group %s', length(nnx), ifelse(xl$label=='', '[ALL]', xl$label))
    hh = pairing[nnx, ]

    myCor = computePairedCorrelation(matrixA=aMtx[, xl$ccn], matrixB=bMtx[, xl$ccn], pairing=hh, BPP=BPP)
    set(myCor, j=c('group', 'aIdx', 'bIdx'), value=list(group=rep(ifelse(xl$label=='', '[ALL]', xl$label), nrow(myCor)), aIdx=unlist(hh[, 1]), bIdx=unlist(hh[, 2])))
    myCorrel = rbind(myCorrel, myCor)
    
    if(doBiweight) {
      aBcw = Matrix::t(Matrix::Matrix(WGCNA::bicovWeights(t(aMtx[, xl$ccn])), sparse=T))
      bBcw = Matrix::t(Matrix::Matrix(WGCNA::bicovWeights(t(bMtx[, xl$ccn])), sparse=T))
      myCorB = computePairedCorrelation(matrixA=aMtx[, xl$ccn], matrixB=bMtx[, xl$ccn], pairing=hh, weightA=aBcw, weightB=bBcw, BPP=BPP)
      set(myCorB, j=c('group', 'aIdx', 'bIdx'), 
          value=list(group=paste0('bw.', rep(ifelse(xl$label=='', '[ALL]', xl$label), nrow(myCor))), aIdx=unlist(hh[, 1]), bIdx=unlist(hh[, 2])))
      myCorrel = rbind(myCorrel, myCorB)
    }
    flog.info('Completed correlations on pairs in group %s', ifelse(xl$label=='', '[ALL]', xl$label))
  }
  return(myCorrel)
}

# Compute paired correlations driver functions

# Get strand aware combination pairs from two ranges like annotations, within a distance.
# If stranded = T, a column strand is expected in the B annotation frame, and pairs will only be computed based on the TSS site 
# (TSS=start if strand='+' or TSS=end if strand='-')
# If stranded = F pairs will be computed to any overlap between intervals of A and B up to distance to each side.
# 
# A and B are expected to have a column 'idx' with unique values. The function returns a data table with aIdx and bIdx pairs
getFeaturePairs <- function(annoA, annoB, stranded=TRUE, distance=1000) {
  
  cna = c('seqnames', 'start', 'end', 'idx')
  stopifnot(all(cna %in% colnames(annoA)) && all(cna %in% colnames(annoB)) && (!stranded || 'strand' %in% colnames(annoB)))
  
  # Create extended windows on annoA, up to distance to each side
  mvExt = annoA[, .(seqnames, start=start - distance, end=end + distance, idx)]
  # Check for negative starts
  mvExt[start < 1, start := 1]
  setkey(mvExt, seqnames, start, end)
  
  if (stranded) {
    # Forward stranded TSSs
    txFwdTSS = annoB[strand != '-', .(seqnames, start, end=start, idx)]
    setkey(txFwdTSS, seqnames, start, end)
    ggF = foverlaps(mvExt, txFwdTSS, nomatch=NULL, which=TRUE)
    ggF = data.table(aIdx = mvExt[['idx']][ggF[['xid']]],
                     bIdx = txFwdTSS[['idx']][ggF[['yid']]])
    
    # Reverse stranded TSSs  
    txRevTSS = annoB[strand != '+', .(seqnames, start=end, end, idx)]
    setkey(txRevTSS, seqnames, start, end)
    ggR = foverlaps(mvExt, txRevTSS, nomatch=NULL, which=TRUE)
    ggR = data.table(aIdx = mvExt[['idx']][ggR[['xid']]],
                     bIdx = txRevTSS[['idx']][ggR[['yid']]])
    
    opairs = unique(rbind(ggF, ggR))
    
  } else {
    
    txG = annoB[, .(seqnames, start, end, idx)]
    setkey(txG, seqnames, start, end)
    ggX = foverlaps(mvExt, txG, nomatch=NULL, which=TRUE)
    
    opairs = data.table(aIdx = mvExt[['idx']][ggX[['xid']]],
                        bIdx = txG[['idx']][ggX[['yid']]])
  }
  setkey(opairs, aIdx, bIdx)
  return(opairs)
}

