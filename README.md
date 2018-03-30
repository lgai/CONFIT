# CONFIT - association testing in multiple traits

The code for CONFIT has been split into 4 parts:
  1. simulate or read in GWAS summary statistics, and compute CONFIT test statistic.
  
  2. get null distribution of CONFIT test statistic. (these should be run in parallel)
  
  3. get top null CONFIT test statistics (used to get high resolution p-values for significant/near significant SNPs).
  
  4. get p-values for CONFIT

