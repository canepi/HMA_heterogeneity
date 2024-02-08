HMA Heterogeneity
=========
This repo contains source code for the manuscript ['Upregulated cholesterol biosynthesis facilitates the survival of methylation-retaining AML cells following decitabine treatment'.](https://www.biorxiv.org/content/10.1101/2024.01.30.577864v1)


Abstract
--------
DNA hypomethylating agents (HMAs) are used to treat acute myeloid leukaemia (AML) and myelodysplasia patients who are unsuitable for intensive chemotherapy, but low response rates and therapy-resistant relapse remain significant challenges. To optimise HMA efficacy, we must understand how resistance and relapse arise from cells that survive treatment. Here we combine single-cell multi-omic analysis with parallel colony-forming assays to link HMA-induced molecular heterogeneity with functional consequences in AML cells. HMAs, azacytidine (AZA) and decitabine (DAC), induced global epigenetic heterogeneity associated with upregulation of inflammatory responses and cell death pathways in a subset of hypomethylated cells. Some AML cells maintained high DNA methylation during treatment, and these methylation-retaining cells had increased self-renewal capacity following DAC, but not AZA. Molecular profiling of individual colonies revealed upregulated cholesterol biosynthesis as an adaptation to HMA treatment, and inhibition by rosuvastatin enhanced DAC effects in vitro and in vivo. Thus, HMA-induced heterogeneity has important implications for AML cell growth and statins are a candidate co-treatment strategy to delay or prevent HMA-resistant relapse.

<p align="center">
 <img src="../230815_Figure1.png" style="width: 30%; height: 30%"/>​
</p>



Content:
--------
This directory contains the source code for the prepocessing, normalisation and batch correction of NOMe-seq and RNA-seq data

* `/scNOMe/`: Pre-processing, normalisation, batch correction and rate calculation of NOMe-seq data.
* `/scRNA/`: Pre-processing, Normalisation and batch correction of RNA-seq count data.


Data:
--------
The raw data is accessible at GEO [GSEXXXXXX](add_hyperlink).  


Acknowledgements/Contact:
--------
Manuscript correspondence: Dr Danielle Bond (danielle.bond@newcastle.edu.au) and Dr Heather Lee (heather.lee@newcastle.edu.au)  
Computational analysis: Dr Sean Burnard (sean.burnard@newcastle.edu.au)

# To Do:
1. Add GEO accession ID and hyperlink after upload completion and release.