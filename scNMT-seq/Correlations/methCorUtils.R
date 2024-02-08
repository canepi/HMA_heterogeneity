#
# More utility functions taken from package I wrote bu now defunct
#
# CR 20240207
##
#
# Takes a regions definition file and returns data.table, sorted
get_regions_file(
  aFile,
  colInfo = c(chr = 1, start = 2, end = 3),
  genome = "GRCh38",
  zero_based = c("none", "start", "end", "both"),
  minWinLen = NULL,
  maxWinLen = NULL
)
{
  stopifnot(all(c("chr", "start", "end") %in% names(colInfo)))
  zero_based = zero_based[1]
  dd = data.table::fread(file = aFile)
  oldcn = colnames(dd)[colInfo]
  data.table::setnames(dd, oldcn, names(colInfo))
  if (zero_based %in% c("start", "both")) {
    dd[, `:=`(start, start + 1)]
  }
  if (zero_based %in% c("end", "both")) {
    dd[, `:=`(end, end + 1)]
  }
  dd[, `:=`(chr, gsub("^chr", "", as.character(chr)))]
  data.table::setkey(dd, chr, start, end)
  return(dd)
}


# Get a standard parameters for use in methylation correlation scripts
methRegionsPar <- function (...)
{
  rl <- list(
    parallel = TRUE,
    numWorkers = 5,
    parallelLog = TRUE,
    logThreshold = "INFO",
    cellWeights = "var",
    logDir = NULL,
    jobname = "METHCOR"
  )
  aa <- list(...)
  for (nn in names(aa)) {
    if (nn %in% names(rl)) {
      rl[[nn]] <- aa[[nn]]
    }
  }
  return(rl)
}

#
# Get a BiocParallel object
getBPP <- function (params = methRegionsPar())
{
  bpp <- NULL
  if (params[["parallel"]]) {
    bpp <-
      BiocParallel::MulticoreParam(workers = params[["numWorkers"]],
                                   log = params[["parallelLog"]],
                                   threshold = params[["logThreshold"]])
  }
  else {
    bpp <- BiocParallel::SerialParam(log = params[["parallelLog"]],
                                     threshold = params[["logThreshold"]])
  }
  if (!is.null(params[["logDir"]])) {
    BiocParallel::bplogdir(bpp) <- params[["logDir"]]
  }
  if (!is.null(params[["jobname"]])) {
    BiocParallel::bpjobname(bpp) <- params[["jobname"]]
  }
  else {
    BiocParallel::bpjobname(bpp) <- "SAMPM"
  }
  return(bpp)
}

# Get a standard regions genome-wide, with length winLen and step winStep
get_all_std_chr_regions <- function (winLen = NULL, winStep = NULL)
{
  dd = data.table::fread(
    text = "
      chr  start end bases ID strand
      1     1   250000000   248,387,328   chr1   *
      2     1   245000000   242,696,752   chr2   *
      3     1   202000000   201,105,948   chr3   *
      4     1   200000000   193,574,945   chr4   *
      5     1   190000000   182,045,439   chr5   *
      6     1   180000000   172,126,628   chr6   *
      7     1   170000000   160,567,428   chr7   *
      8     1   150000000   146,259,331   chr8   *
      9     1   160000000   150,617,247   chr9   *
      10    1   140000000   134,758,134   chr10  *
      11    1   140000000   135,127,769   chr11  *
      12    1   135000000   133,324,548   chr12  *
      13    1   115000000   113,566,686   chr13  *
      14    1   105000000   101,161,492   chr14  *
      15    1   101000000    99,753,195   chr15  *
      16    1   100000000    96,330,374   chr16  *
      17    1    88000000    84,276,897   chr17  *
      18    1    85000000    80,542,538   chr18  *
      19    1    65000000    61,707,364   chr19  *
      20    1    69000000    66,210,255   chr20  *
      21    1    50000000    45,090,682   chr21  *
      22    1    55000000    51,324,926   chr22  *
      X     1   160000000   154,259,566   chrX   *
      Y     1    65000000    62,460,029   chrY   *
      MT    1       20000        16,569   chrM   *"
  )
  
  dd = data.table::setkeyv(dd, cols = c("chr", "start", "end"))
  if (!is.null(winLen)) {
    winStep = if (is.null(winStep) || winStep <= 0) winLen else winStep
    ll = NULL
    for (j in seq(nrow(dd))) {
      ll = c(ll,
             lapply(seq(
               from = dd[j, start],
               to = dd[j, end] + winStep - 1,
               by = winStep
             ),
             function(stp, wLen, achr, aid) {
               list(
                 chr = achr,
                 start = stp,
                 end = stp + wLen - 1,
                 bases = wLen,
                 ID = paste0(aid, ".", stp),
                 strand = "*"
               )
             }, achr = dd[["chr"]][j], wLen = winLen, aid = dd[["ID"]][j]))
    }
    dd2 = data.table::rbindlist(ll)
    dd = dd2
  }
  data.table::set(dd, j = "UID", value = list(dd[, sprintf("%s:%d-%d", chr, start, end)]))
  return(dd)
}

methVar_dataset <- function (
    sampleFiles,
    methWindows,
    export_all = FALSE,
    sampleNames = sapply(sampleFiles, function(x) unlist(strsplit(basename(x[1]), ".", fixed = T))[1], USE.NAMES = F),
    alpha = 1,
    beta = 1,
    params = methRegionsPar()
)
{
  stopifnot(length(sampleFiles) == length(sampleNames))
  stopifnot(inherits(methWindows, "data.table"))
  if (is.null(params))
    params = methRegionsPar()
  export_cols = c("rate", "rate.var", "loci.n")
  if (export_all) {
    export_cols = c(export_cols, "meth.n", "unmeth.n", "total.n")
  }
  flog.info("Begin process of samples methylation data")
  bpp = getBPP(params = params)
  mRegions = methWindows
  mRegions[, `:=`(mvdsID, .I)]
  res <- BiocParallel::bplapply(
    sampleFiles,
    methVar_regions,
    summaryIntervals = mRegions,
    alpha = alpha,
    beta = beta,
    BPPARAM = bpp
  )
  flog.info("Finished per window and sample process (parallel: %s)",
            params[["parallel"]])
  flog.info("Merging results")
  for (k in seq_along(res)) {
    mRegions = cbind(mRegions, res[[k]]$data[mRegions, .SD,
                                             .SDcols = export_cols, on = "mvdsID"])
    data.table::setnames(mRegions,
                         old = export_cols,
                         new = paste(sampleNames[k],
                                     export_cols, sep = "."))
  }
  flog.info("Done merging")
  mRegions[, `:=`(mvdsID, NULL)]
  mRegions = mRegions[chr %in% c(as.character(seq(22)), "X",
                                 "Y", "MT"),]
  return(mRegions)
}


methVar_regions <- function (fileList, summaryIntervals, alpha = 1, beta = 1) 
{
  summaryIntervals = data.table::setDT(summaryIntervals)
  oldkey = data.table::key(summaryIntervals)
  data.table::setkey(summaryIntervals, chr, start, end)
  z1 = microbenchmark::microbenchmark({
    ll <- list()
    ll = lapply(fileList, function(aFile, smry) {
      futile.logger::flog.debug("Process coverage file %s", 
                                aFile)
      tt = data.table::fread(file = aFile)
      data.table::setnames(tt, c("chr", "start", "end", 
                                 "meth.pc", "meth.n", "unmeth.n"))
      tt[, `:=`(chr, as.character(chr))]
      data.table::setkey(tt, chr, start, end)
      hh = data.table::foverlaps(smry, tt, nomatch = NULL, 
                                 which = T)
      methsum = tt[hh$yid, ][, `:=`("_iv", hh$xid)]
      methsum = methsum[, .(meth.n = sum(meth.n, na.rm = T), 
                            unmeth.n = sum(unmeth.n, na.rm = T), loci.n = .N), 
                        by = "_iv"]
      methsum
    }, smry = summaryIntervals)
    gp = data.table::rbindlist(ll)
    gp = gp[, lapply(.SD, sum), by = "_iv"]
    gp = cbind(summaryIntervals[gp$"_iv", ], gp)
    data.table::set(gp, j = c("_iv", "total.n"), value = list(NULL, 
                                                              gp[, meth.n] + gp[, unmeth.n]))
    data.table::setkeyv(gp, cols = oldkey)
  }, times = 1)
  rm(ll)
  gc(verbose = F)
  futile.logger::flog.info("Files in:\n%s", paste0(fileList, 
                                                   collapse = "\n"))
  futile.logger::flog.info("Regions: %d, time: %ss", nrow(gp), 
                           z1$time * 1e-09)
  gp[, `:=`(rate, (meth.n + alpha)/(total.n + alpha + beta))]
  gp[, `:=`(rate.var, rate * (1 - rate)/total.n)]
  return(list(data = gp, files = fileList, time = z1$time))
}


SumExpFromRateFile <- function (rateFile = NULL, methData = NULL, samples = NULL, colData = NULL, excludeRowInfo = "", joinchar = ".") 
{
  if (!xor(is.null(rateFile), is.null(methData))) {
    stop("Only one of 'rateFile' or 'methData' must be given")
  }
  else if (!is.null(rateFile)) {
    require(data.table)
    methData = fread(file = rateFile)
  }
  if (is.null(samples) && !is.null(colData)) {
    samples = rownames(colData)
  }
  if (is.null(samples)) {
    nn = grep(".*[_.]loci.n$", colnames(methData))
    samples = unique(gsub("^(.*)[_.]loci.n$", "\\1", colnames(methData)[nn]))
  }
  if (!is.null(colData) && length(samples) != nrow(colData)) {
    stop("Samples length and column data rows must be same number")
  }
  rdc = unique(unlist(sapply(samples, function(x) grep(x, colnames(methData)), 
                             USE.NAMES = F)))
  rdc = setdiff(seq(ncol(methData)), rdc)
  if (!is.null(colData)) {
    colData = as.data.table(colData)
    set(colData, j = c("idx", "sample"), value = list(idx = seq(nrow(colData)), 
                                                      sample = samples))
  }
  else {
    colData = data.table(sample = samples, idx = seq(length(samples)))
  }
  dt = methData[, .SD, .SDcols = paste0(samples, excludeRowInfo, 
                                        joinchar, .coltypes[["loci"]])]
  irx = which(dt[, Reduce("+", lapply(.SD, is.na))] == ncol(dt))
  if (length(irx) > 0) {
    methData = methData[-irx, ]
  }
  myAssays = list()
  for (asy in names(.coltypes)) {
    dt = methData[, .SD, .SDcols = paste0(samples, excludeRowInfo, 
                                          joinchar, .coltypes[asy])]
    colnames(dt) = samples
    myAssays = c(myAssays, list(dt))
  }
  myAssays = S4Vectors::SimpleList(myAssays)
  names(myAssays) = names(.coltypes)
  if (!is.null(excludeRowInfo) && nchar(excludeRowInfo) > 0) {
    xcl = grep(excludeRowInfo, colnames(methData)[rdc])
    if (length(xcl) > 0) 
      rdc = rdc[-xcl]
  }
  rdat = methData[, .SD, .SDcols = rdc]
  if (!all(rdat$strand %in% c("+", "-", ".", "*"))) {
    rdat[!strand %in% c("+", "-", ".", "*"), `:=`(strand, 
                                                  ".")]
  }
  rdat = GenomicRanges::makeGRangesFromDataFrame(df = rdat, 
                                                 keep.extra.columns = T)
  exps = SummarizedExperiment::SummarizedExperiment(assays = myAssays, 
                                                    rowRanges = rdat, colData = colData)
  return(exps)
}

