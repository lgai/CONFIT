## CONFIT - association testing in multiple traits

The code for CONFIT has been split into 4 parts:
1. simulate or read in GWAS summary statistics, and compute CONFIT test statistic.
  
2. get null distribution of CONFIT test statistic. (these should be run in parallel)
  
3. get top null CONFIT test statistics (used to get high resolution p-values for significant/near significant SNPs).
  
4. get p-values for CONFIT


## Example commands
Below, we show some example commands with options to be used with a small example dataset. 

#### preliminaries:
- requires Python 2.7
- requires numpy and scipy


#### Step 1: Compute the CONFIT test stat on your chosen data

For pylmm-formatted data
python src/confit1_data.py --exptName NFBCexample --outdir exampleOutput  --nSnp 100  --gwasFormat pylmm  --traits NFBC_hdlres.pylmm.out NFBC_ldlres.pylmm.out NFBC_tgres.pylmm.out   --traitPath your/path/to/data/CONFIT_examples/gwasExampleData

For UKBB data in Neale group format
- this sample data is already sorted so that SNPs are in same order for each trait
```
python src/confit1_data.py --exptName UKBBexample --outdir exampleOutput --nSnp 100 --gwasFormat UKBB_sorted  --traits UKBB_6177_1.assoc.tsv.sorted UKBB_6177_2.assoc.tsv.sorted UKBB_6177_3.assoc.tsv.sorted   --traitPath your/path/to/data/CONFIT_examples/gwasExampleData 
```


For simulated data (100 snps, 3 traits, with default settings)
```
python src/confit1_data.py --exptName simulatedExample --outdir exampleOutput --useSimulatedData 1 --nSnp 100 --nTraits 3 --sigmasq_mu_sim 25
```

For simulated data (1000 snps, 3 traits, with optional arguments shown)
```
python src/confit1_data.py --exptName simulatedExampleWithOptions --outdir exampleOutput --useSimulatedData 1 --nSnp 1000 --nTraits 3 --sigmasq_mu_sim 25 --t1trueSigmasq 25  --Sigma_e_file_sim your/path/to/data/CONFIT_examples/simulationExampleData/threeTraitExample_Sigma_z.txt --truePriorFile your/path/to/data/CONFIT_examples/simulationExampleData/threeTraitExample_truePriorFile.txt
```


#### Step 2: Null simulations 
- Note: we'll use the UKBB example data for the rest of the example steps
- In this example, we're running the null simulation step twice, sequentially. In practice, you probably want to do many of these in parallel.
- We do not compute p-values with much precision in the examples, since one would need to generate many null panels to do so.

```
python src/confit2_nullsim4pval.py --exptName UKBBexample --outdir exampleOutput --nulldir exampleNullsim --nTraits 3  --nNullPerRound 4000 --nNullRoundsPerJob 25 --taskID 1

python src/confit2_nullsim4pval.py --exptName UKBBexample --outdir exampleOutput --nulldir exampleNullsim --nTraits 3  --nNullPerRound 4000 --nNullRoundsPerJob 25 --taskID 2
```


#### Step 3: Get top test statistics (which will be used to set the p-value threshold)
```
python src/confit3_getTopStatistics.py --exptName UKBBexample --nulldir exampleNullsim --outdir exampleOutput --nTopNull 1000 --nNullFiles 2 --nNullTotal 20000
```


#### Step 4: Compute p-values for MI GWAS and CONFIT test statistics
```
python src/confit4_consolidatePvals.py --exptName UKBBexample --outdir exampleOutput --nulldir exampleNullsim --useSimulatedData 0 --nTraits 3 --nNullLowRes 5000
```






