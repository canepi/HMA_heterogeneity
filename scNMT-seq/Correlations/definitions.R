#
# Definition of common paths, filenames and parameters used on various scripts and mardown 
# documents in this folder
#
# CR 20201029, 20210121, 20220825
#
projDir = path.expand('/home/carlos/Projects/BIO19003-201707-Canepi-Lee/')
dataDir = paste0(projDir, 'data/20220823-CorrelationAnalysis/')
qcDir   = paste0(projDir, 'data/20220928-QC-3Batches/')
workDir = paste0(projDir, 'docs/20221031-CorrelationAnalysis-3Batches-v2/')
workDir2= paste0(projDir, 'docs/20220909-PWMethDistAnalysis-3Batches/')
workDir3= paste0(projDir, 'docs/20220825-CorrelationAnalysis-3Batches/')
annoDir = paste0(dataDir, 'Annotations/')
genoDir = '/home/bioinf/genomes/hs/ensemble/'
covDir  = paste0(dataDir, 'CovFiles/')
# rnaDir  = paste0(dataDir, 'RNASeq/')
rnaDir  = paste0(projDir, 'data/20230328-RNASeq-3Batches/')

genoSmryFile = paste0(genoDir, 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.summary')

myContext = "GpC"
doBiweight = FALSE

# Use the names of the list to select type
# Is important the filtering fields are uniform names in all datasets:
# The colInfo field indicates which columns in the original annotation file correspond to the key fields
# The names 'chr', 'start', 'end', 'strand', 'ID' are mandatory
# adjustStart is TRUE if the start value is zero-based
# Adding the maximum number of rows to use on partial distance computations, as this is somewhat dataset-dependant
# Must be a valid integer for display purposes. Enough to be larger than the number of regions to do a single step.
featureData = list(
  HL60.H3K4me3 = list(dataFile = paste0(annoDir, 'HL60_H3K4me3_ENCFF021JBH.txt'), 
                      adjustStart = F,
                      colInfo  = c(chr=1, start=2, end=3, ID=5, strand=4),
                      distance = 10000, 
                      match = 'byPos', 
                      step.max.regions = 1000000L),
  HL60.H3K27ac = list(dataFile = paste0(annoDir, 'HL60_H3K27ac_ENCFF763UAG.txt'), 
                      adjustStart = F,
                      colInfo  = c(chr=1, start=2, end=3, ID=5, strand=4),
                      distance = 10000, 
                      match = 'byPos', 
                      step.max.regions = 1000000L),
  CGI_SeqMonk  = list(dataFile = paste0(annoDir, 'CGI_SeqMonk.txt'),  
                      adjustStart = F,
                      colInfo  = c(chr=1, start=2, end=3, ID=5, strand=4),
                      distance = 10000, 
                      match = 'byPos', 
                      step.max.regions = 1000000L),
  Promoter     = list(dataFile = paste0(annoDir, 'Promoter_mRNA_upstream_1500_500.txt'),  
                      adjustStart = F,
                      colInfo  = c(chr=1, start=2, end=3, ID=5, strand=4),
                      distance = 10000, 
                      match = 'byID', 
                      step.max.regions = 1000000L),
  Promoter50   = list(dataFile = paste0(annoDir, 'Promoter_mRNA_upstream_50_50.txt'),  
                      adjustStart = F,
                      colInfo  = c(chr=1, start=2, end=3, ID=5, strand=4),
                      distance = 10000, 
                      match = 'byID', 
                      step.max.regions = 1000000L),
  Promoter4k   = list(dataFile = paste0(annoDir, 'Promoter_mRNA_upstream_2000_2000.txt'),  
                      adjustStart = F,
                      colInfo  = c(chr=1, start=2, end=3, ID=5, strand=4),
                      distance = 10000, 
                      match = 'byID', 
                      step.max.regions = 1000000L),
  Intergenic   = list(dataFile = paste0(annoDir, 'Intergenic2.txt'), 
                      adjustStart = F,
                      colInfo  = c(chr=1, start=2, end=3, ID=5, strand=4),
                      distance = 10000, 
                      match = 'byPos', 
                      step.max.regions = 1000000L),
  Introns      = list(dataFile = paste0(annoDir, 'Introns.txt'), 
                      adjustStart = F,
                      colInfo  = c(chr=1, start=2, end=3, ID=5, strand=4),
                      distance = 10000, 
                      match = 'byID', 
                      step.max.regions = 115000L),
  Exons        = list(dataFile = paste0(annoDir, 'Exons.txt'), 
                      adjustStart = F,
                      colInfo  = c(chr=1, start=2, end=3, ID=5, strand=4),
                      distance = 10000, 
                      match = 'byID', 
                      step.max.regions = 175000L),
  Window3k1k   = list(dataFile = paste0(annoDir, '3kb_window_1kb_step.txt'), 
                      adjustStart = F,
                      colInfo  = c(chr=1, start=2, end=3, ID=5, strand=4),
                      distance = 10000, 
                      match = 'byPos', 
                      step.max.regions = 1550000L)
)
