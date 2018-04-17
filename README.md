# CONFIT - association testing in multiple traits

The code for CONFIT has been split into 4 parts:
  1. simulate or read in GWAS summary statistics, and compute CONFIT test statistic.
  
  2. get null distribution of CONFIT test statistic. (these should be run in parallel)
  
  3. get top null CONFIT test statistics (used to get high resolution p-values for significant/near significant SNPs).
  
  4. get p-values for CONFIT


Below, we show some example commands with options to be used with our small example dataset. (We do not compute p-values with much precision in the examples, since one would need to generate many nul panels to do so.)

############################################################
# preliminaries:
# - make sure you are using Python 2.7
# - make sure you have numpy and scipy

############################################################
# step 1: Compute the CONFIT test stat on your chosen data

#    for pylmm-formatted data.
#    Note: for real data, user must specify the summary statistic format as either "pylmm" or "UKBB_sorted", and number of SNPs in the dataset
python code/confit1_data.py --exptName NFBCexample --outdir exampleOutput  --nSnp 100  --gwasFormat pylmm  --traits NFBC_hdlres.pylmm.out NFBC_ldlres.pylmm.out NFBC_tgres.pylmm.out   --traitPath /your/path/to/gwasExampleData

#    for UKBB data in Neale group format. 
#    Note: this sample data is already sorted so that SNPs are in same order for each trait
python code/confit1_data.py --exptName UKBBexample --outdir exampleOutput --nSnp 100 --gwasFormat UKBB_sorted  --traits UKBB_6177_1.assoc.tsv.sorted UKBB_6177_2.assoc.tsv.sorted UKBB_6177_3.assoc.tsv.sorted   --traitPath /your/path/to/gwasExampleData  --Sigma_Z_file /your/path/to/exampleOutput/UKBB_Sigma_Z.txt

#    for simulated data (100 snps, 3 traits)
python code/confit1_data.py --exptName simulatedExample --outdir exampleOutput --useSimulatedData 1 --nSnp 100 --nTraits 3 --sigmasq_mu_sim 25

############################################################
# step 2: Null simulations 
# (Note: we'll use the UKBB example data for the rest of the example steps)
# (Note: In these example, we're running the null simulation step twice, sequentially, with a different taskID for each output file. In real life, you probably want to do many of these in parallel)

python code/confit2_nullsim4pval.py --exptName UKBBexample --nulldir exampleNullsim --nTraits 3   --Sigma_Z_file /your/path/to/exampleOutput/UKBB_Sigma_Z.txt    --wtsFile /your/path/to/exampleOutput/UKBBexample_wt.txt --nNullPerRound 4000 --nNullRoundsPerJob 25 --taskID 1

python code/confit2_nullsim4pval.py --exptName UKBBexample --nulldir exampleNullsim --nTraits 3   --Sigma_Z_file /your/path/to/exampleOutput/UKBB_Sigma_Z.txt    --wtsFile /your/path/to/exampleOutput/UKBBexample_wt.txt --nNullPerRound 4000 --nNullRoundsPerJob 25 --taskID 2

############################################################
#    step 3: Get top test statistics (which will be used to set the p-value threshold)

python code/confit3_getTopStatistics.py --exptName UKBBexample --nulldir exampleNullsim --outdir exampleOutput --nTopNull 1000 --nNullFiles 2 --nNullTotal 20000

############################################################
# step 4: Compute p-values for MI GWAS and CONFIT test statistics

python code/confit4_consolidatePvals.py --exptName UKBBexample --outdir exampleOutput --nulldir exampleNullsim --useSimulatedData 0 --nTraits 3 --nNullLowRes 5000





