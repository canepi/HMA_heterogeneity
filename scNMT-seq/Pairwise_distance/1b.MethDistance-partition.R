# Pairwise Cell Methylation Distance analysis over regions
#
# 1. Computation of region-pair distances and summary over regions
#
# Regions as defined for the analysis of 3Batches methylation  data
#
# CR 220909 based on preliminary analysis on May/June 2022
#
# Some code is shared with Correlation analysis
#
# This implementation attempts to use Apache/Arrow datasets
library(data.table)
library(arrow)
library(dplyr)
library(futile.logger)
# if(!require(R.utils)) {
#   install.packages('R.utils')
# }

if(!require(methCor)) {
  install.packages('/home/carlos/Workspace/methCor_0.0.1.1_R_x86_64-redhat-linux-gnu.tar.gz', repos=NULL)
  Sys.sleep(3)
  library(methCor)
}

### Set up ####
source('/home/carlos/Projects/BIO19003-201707-Canepi-Lee/docs/20220909-PWMethDistAnalysis-3Batches/definitions.R')
source('/home/carlos/Projects/BIO19003-201707-Canepi-Lee/docs/20220909-PWMethDistAnalysis-3Batches/methdist-utils.R')

# Could not find this published or a way to deserialise from a dataset. See https://rdrr.io/cran/arrow/src/R/metadata.R
.deserialise_arrow_r_metadata <- function(x) {
  tryCatch(expr = {
    out <- unserialize(charToRaw(x))
    # if this is still raw, try decompressing
    if (is.raw(out)) {
      out <- unserialize(memDecompress(out, type = "gzip"))
    }
    out
  }, error = function(e) { 
    warning("Invalid metadata$r", call. = FALSE)
    NULL
  })
}

do.DEBUG = TRUE
do.forceGLB = FALSE   # Force computation of per cell globals

outDir = paste0(projDir, 'docs/20220910-PWMethDist-3B/')
# This is the output directory where analysis will be saved
if(!dir.exists(outDir)) {
  dir.create(outDir, showWarnings=F, recursive=T)
}
setwd(outDir)

logFile = strftime(Sys.Date(), "methdist-%Y%m%d.log")
flog.appender(flog.indent_wrapper(appender=appender.tee(logFile), indent="-       "))
if(do.DEBUG)
  flog.threshold(DEBUG)

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

flog.info("####\n\t\t   #### START Region methylation distance calculation\n\t\t   ####")

mvpar = methCor::methRegionsPar(logDir=outDir, numWorkers=10, jobname='Dist')

BPP = getBPP(mvpar)

featureData = loadAnnotationData(featureData, add.byChr=F, BPP=BPP)

#### Methylation data ######

# This is the main processing loop #####################
myContext = 'GpC'
for(myContext in c('GpC','CpG')) {
  
  flog.info("----\n\t\t   ---- START Process for context %s\n\t\t   ----", myContext)
  
  covFiles = list.files(path=covDir, pattern=paste0(".*.NOMe.", myContext, ".cov.gz"), full.names=T, recursive=F)
  cvnames = .inferColData(basename(covFiles))[, paste0(Treatment,'.',Batch,'.',Cell)]
  names(covFiles) = cvnames
  oo = order(cvnames)
  covFiles = covFiles[oo]
  cvnames = cvnames[oo]
  
  globalSmry = NULL
  globalCellsInfo = NULL
  featName = "Window3k1k"
  # "HL60.H3K4me3" "HL60.H3K27ac" "CGI_SeqMonk" "Promoter" "Promoter50" "Promoter4k" "Intergenic" "Introns" "Exons" "Window3k1k" "All_loci" 
  for(featName in names(featureData)) {
    FD = featureData[[featName]]
    flog.info('>>>>>> Start process of %s %s annotation', myContext, featName)
    
    cellsInfo = NULL
    ctxSummaryPW = NULL
    gc()

    outfname = paste0(outDir, 'summaryPW-', featName, '-', myContext, '.d')
    ciglb_outfname = gsub('\\.d$', '-cellsInfo.rds', outfname)
    
    ### Global per cell annotation globals
    if(!file.exists(ciglb_outfname) || do.forceGLB) {
      
      cellsInfo = globalMeansByCell(FD[['data']], covFiles, BPPARAM=BPP)
      set(cellsInfo, j=c('Annotation','Context','Step'), value=list(Annotation=featName, Context=myContext, Step='global'))
      saveRDS(cellsInfo, file=ciglb_outfname)
      flog.info('Save by cell Info to %s', ciglb_outfname)
      
    } else {
      
      cellsInfo = readRDS(ciglb_outfname)
      flog.info('Read by cell Info from %s', ciglb_outfname)
      
    }
    globalCellsInfo = rbind(globalCellsInfo, cellsInfo)
    
    ## Pairwise distance computation
    l1names = cvnames[-length(cvnames)]
    partStatus = checkArrowDataset(path=outfname, cvnames=l1names, totalSteps=nrow(FD[['data']]))
    flog.info('Precompute Arrow dataset status:', partStatus, capture=T)
    if(is.null(partStatus) || 
       (nrow(partStatus) > 1 && !all(is.na(partStatus[['L2.complete']])) && nrow(partStatus[(L1.complete) & L2.last > 0, ]) < length(l1names))) {
      # If the dataset is found to have multiple step sizes (difference between consecutive steps, whether for some cells or for all)
      # then NULL is returned. I guess more effort could be put here, but I don't see it being practical.
      flog.warn('Error in dataset %s. SKIP computation', basename(outfname))
      next
    } else if(nrow(partStatus) == 0 ||  # Does not exist
              (nrow(partStatus) == length(l1names) && all(is.na(partStatus[['L2.complete']])) && nrow(partStatus[(L1.complete),]) < length(l1names)) || # Incomplete L1
              (nrow(partStatus) == length(l1names) && !any(is.na(partStatus[['L2.complete']])) && nrow(partStatus[(L2.complete), ]) < length(l1names)) ) {  # Incomplete L2
      # If the dataset is not found at `path` the returned table has 0 rows.
      # If the dataset does not have L1 partition the returned table has 1 row with (L1.name='', L1.complete=T, L2.complete=NA, L2.last=NA)
      # If the dataset has only L1 partition, the returned table has all the L1.names, with L1.complete corresponding to whether the
      # corresponding partition is found.
      # If the dataset is found to have a L2 partition with Step=1, then is assumed to be "old style" where the step number was the first step
      # (not the last) in the partition. In this case, L2.complete is NA and L2.last = 0 as we can not know what the last step was.
      # If the step size is uniform, L2.last will indicate the last step found for the cell and L2.complete whether L2 is done or not.
      toDo = NULL
      if(nrow(partStatus) == length(l1names) && all(is.na(partStatus[['L2.complete']])) && nrow(partStatus[(L1.complete), ]) < length(l1names)) {
        toDo = unlist(partStatus[!(L1.complete), L1.name])
      }

      nstep = ifelse(is.null(FD[['step.max.regions']]), nrow(FD[['data']]), FD[['step.max.regions']])
      stepi = if(nrow(partStatus) > 0) partStatus[, min(L2.last)] + 1 else 1
      while(stepi < nrow(FD[['data']])) {
        stepe = min(stepi + nstep - 1, nrow(FD[['data']]))
        # These the defaults if no steps necessary
        step.part = list()
        step.smry = list(global=NULL, region='UID')
        if(nstep < nrow(FD[['data']])) {
          # These the modified values for the step
          step.part = list('Step' = stepe)
          names(step.smry)[1] = paste0('Step_', stepe)    # We just want the name of the summary for the step
        }
        numRegions = pairwise_dists_4_ctx(FD[['data']][seq(stepi, stepe), ], covFiles, outDir=outfname, do.partitions=toDo, 
                                          other.partitions=step.part, by.smry=step.smry, cov.workers=FD[['cov.workers']], BPPARAM=BPP)

        flog.info('  Completed step %d ', stepe)
        stepi = stepi + nstep
      }
      flog.info('Save Distance to arrow dataset %s', basename(outfname))
    } else {
      flog.info('Using existing Distance arrow dataset %s', basename(outfname))
    }
    
    partStatus = checkArrowDataset(path=outfname, cvnames=l1names, totalSteps=nrow(FD[['data']]))
    flog.info('Postcompute Arrow dataset status:', partStatus, capture=T)
    if(nrow(partStatus) == 1 || 
       (all(is.na(partStatus[['L2.complete']])) && nrow(partStatus[(L1.complete), ]) == length(l1names)) ||
       (nrow(partStatus[(L2.complete), ]) == length(l1names)) ) {
      
      # Get or compute globals per cell pair
      ctxSummaryPW = open_dataset(outfname)
      flog.info('Open Distance arrow file %s', basename(outfname))
      flog.info('', ctxSummaryPW, capture=T)

      # If level 2 dataset, compute from steps
      if(nrow(partStatus[(L2.complete), ]) == length(l1names)) {
        ctxGlobal = ctxSummaryPW %>% 
          filter(grepl('Step_', Summary, fixed=T)) %>% 
          group_by(Cell_A, Cell_B) %>% 
          summarise(pca = sum(`Common cytosines` * `pc meth A`),
                    pcb = sum(`Common cytosines` * `pc meth B`), 
                    n = sum(`Common cytosines`),
                    nr = sum(`Common regions`),
                    mmd = sum(`Common cytosines` * `Mean meth distance`),
                    dmm = mean(`Delta mean meth`)
          ) %>% 
          mutate(`Summary` = 'global',
                 `pc meth A` = pca / n,
                 `pc meth B` = pcb / n,
                 `Common cytosines` = n,
                 `Common regions` = nr,
                 `Mean meth distance` = mmd / n,
                 `Delta mean meth` = dmm,
                 `Annotation` = featName,
                 pca = NULL,
                 pcb = NULL,
                 n = NULL,
                 nr = NULL,
                 mmd = NULL,
                 dmm = NULL
          ) %>%
          arrange(Cell_A, Cell_B) %>%
          collect()
        
        flog.info('Computing pairwise globals from L2 dataset', ctxGlobal, capture=T)
        globalSmry = rbind(globalSmry, as.data.table(ctxGlobal), fill=T)
        
      } else {
        ctxGlobal = ctxSummaryPW %>% 
          filter(Summary == 'global') %>% 
          select(!UID) %>% 
          mutate(Annotation = featName) %>% 
          collect()
        flog.info('Retrieve pairwise globals from L1 or L0 dataset', ctxGlobal, capture=T)
        globalSmry = rbind(globalSmry, as.data.table(ctxGlobal), fill=T)
      }
      
    } else {
      flog.fatal('Error in dataset %s. STOP', basename(outfname))
      stop(0)
    } 

    flog.info('End distance calculation for %s %s', myContext, featName)
  }
  
  globalSmryFile = paste0(outDir, 'summaryPW-global-', myContext, '.rds')
  saveRDS(as.data.table(globalSmry), file=globalSmryFile)
  flog.info('Saved global summary PW information to %s', globalSmryFile)
  
  globalSmryFile = paste0(outDir, 'cellPW-global-', myContext, '.csv')
  fwrite(globalCellsInfo, file=globalSmryFile)
  flog.info('Saved global cell PW information to %s', globalSmryFile)
  
  rm(globalSmry, globalCellsInfo)
  gc()
  
}
flog.info('%%%%%%% END PROCESS %%%%%%%%')

stop(0)
