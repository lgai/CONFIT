
## CONFIT - association testing in multiple traits

The code for CONFIT has been split into 4 parts:

1. Simulate or read in GWAS summary statistics, and compute CONFIT test statistic for each SNP. 

2. Get the null distribution of CONFIT test statistic. This should be run in parallel, since many null test statistics are needed to get the distribution.

3. Get top null CONFIT test statistics, which are used to help compute p-values.
  
4. Compute CONFIT p-values for each SNP.



## Tutorial
We describe basic use of CONFIT and show sample commands that can be run on the provided example data sets. 

### Preliminaries
- requires Python 2.7
- requires numpy and scipy

### Step 1: Compute the CONFIT test statistics on chosen data
In this step, we read in or simulate GWAS summary statistics, and compute CONFIT test statistic for each SNP from the GWAS data. 

General options 
`exptName` - name for output files (without extension) 
`outdir` - output directory 


#### If using GWAS data from file
Required options for GWAS data from file
 - `nSnp` - number of SNPs in analysis
 - `gwasFormat` - specify format of GWAS summary statistics
- `traits` and `traitPath` - names of each summary statistic file and their directory

Example for pylmm-formatted data
```
python src/confit1_data.py --exptName NFBCexample --outdir exampleOutput  --nSnp 100  --gwasFormat pylmm  --traits NFBC_hdlres.pylmm.out NFBC_ldlres.pylmm.out NFBC_tgres.pylmm.out   --traitPath your/path/to/data/CONFIT_examples/gwasExampleData
```

Example for UKBB data in Neale group format. (Note this sample data has already been sorted so that SNPs are in same order for each trait.)
```
python src/confit1_data.py --exptName UKBBexample --outdir exampleOutput --nSnp 100 --gwasFormat UKBB_sorted  --traits UKBB_6177_1.assoc.tsv.sorted UKBB_6177_2.assoc.tsv.sorted UKBB_6177_3.assoc.tsv.sorted   --traitPath your/path/to/data/CONFIT_examples/gwasExampleData 
```

#### If simulating GWAS data
Required options for simulated data
- `useSimulatedData` - Set as 1 use simulated data
- `nSnp` - how many SNPs to simulate
- `sigmasq_mu_sim` - variance of simulated z-scores

Example for simulated data (100 SNPs, 3 traits)
```
python src/confit1_data.py --exptName simulatedExample --outdir exampleOutput --useSimulatedData 1 --nSnp 100 --nTraits 3 --sigmasq_mu_sim 25
```

Example for simulated data (1000 SNPs, 3 traits, with some additional simulation options shown)
```
python src/confit1_data.py --exptName simulatedExampleWithOptions --outdir exampleOutput --useSimulatedData 1 --nSnp 1000 --nTraits 3 --sigmasq_mu_sim 25 --t1trueSigmasq 25  --Sigma_e_file_sim your/path/to/data/CONFIT_examples/simulationExampleData/threeTraitExample_Sigma_z.txt --truePriorFile your/path/to/data/CONFIT_examples/simulationExampleData/threeTraitExample_truePriorFile.txt
```


### Step 2: Null simulations of CONFIT test statistic
This step get draws from the null distribution of the CONFIT test statistic. The null test statistics are then used to obtain p-values in Steps 3 and 4.

Each time confit2_nullsim4pval is run, it will generate null test statistics in a separate file (so this step is easily parallelized in order to generate many null test statistics). 

Required options
- `exptName` - name for output files, as in Step 1
- `outdir` - as in Step 1
- `nulldir` - where to output the null test statistics (which are just intermediate output)
- `nTraits` OR `traits` - either specify how many traits or list them out as in Step 1
- `taskID` - used to give each null test statistic file a unique name. Specify a number for each run of confit2_nullsim4pval.py (i.e. --taskID 2 for the second run).

Other options
- `nNullPerRound` (default $2*10^6$) 
- `nNullRoundsPerJob` (default $25$) - Each time the command is run, it will generate nNullPerRound*nNullRoundsPerJob null test statistics. If running on a system with limited memory, you can decrease `nNullPerRound` so fewer test statistics are generated at a time.

Example
- Note: we'll use the UKBB example data for the rest of the steps
- The below commands run the null simulation step twice, sequentially. This will create two files, each with 4000*25 null test statistics.
- In practice, you probably want to run many of these in parallel, in order to obtain p-values with more precision.
```
python src/confit2_nullsim4pval.py --exptName UKBBexample --outdir exampleOutput --nulldir exampleNullsim --nTraits 3 --nNullPerRound 4000 --nNullRoundsPerJob 25 --taskID 1

python src/confit2_nullsim4pval.py --exptName UKBBexample --outdir exampleOutput --nulldir exampleNullsim --nTraits 3  --nNullPerRound 4000 --nNullRoundsPerJob 25 --taskID 2
```

### Step 3: Get top null test statistics 
This is an intermediate step before the p-values are computed in Step 4. It helps compute high-resolution p-values for the most significant test statistics. The default settings assume $5*10^9$ null test statistics were generated in Step 2.

Options
- `exptName`,  `outdir`, `nulldir` - as in Step 1 and 2
- `nNullFiles` - how many files were generated in Step 2
- `nNullTotal` (default $5 * 10^9$) - how many null test statistics were generated in Step 2 across all files
- `nTopNull` (default $1000$) - how many of top null test statistics to get

Example 
 - with parameters set to match what we did in Step 2.
```
python src/confit3_getTopStatistics.py --exptName UKBBexample --nulldir exampleNullsim --outdir exampleOutput --nNullFiles 2 --nNullTotal 20000
```


### Step 4: Compute p-values for CONFIT test statistics
Compute CONFIT p-values. The default settings assume $5*10^9$ null test statistics were generated in Step 2.
 
Required options
- `exptName`,  `outdir`, `nulldir` - as in Step 1 and 2
- `nTraits` OR `traits` - as in Step 1 and 2
- `useSimulatedData` - 1 if data was simulated by CONFIT, 0 if read from file (used because the output format will have some additional columns if the data was simulated)

By default, CONFIT computes less-significant p-values with lower precision. These additional options may be used to adjust the p-value resolution.
- `nNullFullRes` (default $5 * 10^9$) - how many null test statistics to use for highest resolution. You can set it to `nNullTotal` as in Step 3.
- `nNullLowRes` (default $10^8$) - upper bound on how many null test statistics to use for lower resolution p-values. (For the toy example below, we just set it to a small value.)

Example
```
python src/confit4_consolidatePvals.py --exptName UKBBexample --outdir exampleOutput --nulldir exampleNullsim --useSimulatedData 0 --nTraits 3 --nNullLowRes 5000 --nNullFullRes 20000
```





