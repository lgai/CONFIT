# First step of CONFIT
# (0) read in options
# (1) either generate or read in plink formatted data
# (2) set priors from data and compute CONFIT test statistic for each snp
# (3) write priors + GWAS and CONFIT results to file


from __future__ import print_function

import sys
import getopt
import argparse
import itertools
import collections
import numpy as np
import numpy.random
import scipy.stats
from scipy.stats import multivariate_normal as mvnorm


from confitFunctions import *

################ PARSE OPTIONS ##################
parser = argparse.ArgumentParser(description="CONFIT")

# general options
parser.add_argument("--exptName", type=str, required=True) # what to name output files
parser.add_argument("--outdir", type=str,default=".") # where to put output files
parser.add_argument("--useSimulatedData", type=int, default=0) # whether using simulated or real data
parser.add_argument("--nSnp",  type=int, required=True) # how many snps to simulate, or how many snps in real data
parser.add_argument("--Sigma_Z_file", type=str, default=None) # name of Sigma_ to use for simulation, or where to write Sigma_Z estimated from Finland data. If none and simulating data, use identity

# simulation options
parser.add_argument("--nTraits",  type=int, default=None)
parser.add_argument("--sigmasq_mu_sim", type=float, default=None)
parser.add_argument("--t1trueSigmasq", type=float, default=None) # if t1 should have different variance of lambda than other traits, set it here
parser.add_argument("--truePriorFile", default=None)

# real data options
parser.add_argument("--traits", type=str, nargs='*', default=None) # file name of each trait, including filetype extension
parser.add_argument("--traitPath",type=str,default=None) # if using real data, must specify directory where gwas summary stats are located
parser.add_argument("--gwasFormat", type=str, default=None) # either use "pylmm" (https://github.com/nickFurlotte/pylmm)
# or use use "UKBB_sorted" (https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/README.md#summary-stat-output), but snps must be sorted to all be in same order

# options for CONFIT test statistic
parser.add_argument("--initialPriorFile", type=str, default=None) # user may specify their own priors on each configuration
parser.add_argument("--sigmasq_mu_guess", type=float, default=25)
parser.add_argument("--avg1active",type=int,default=1) # whether to average over configs w 1 active trait
parser.add_argument("--countPriorThresh",type=float,default=1e-4) # threshold for count prior
parser.add_argument("--minPriorC", type=float, default=1e-20) # minimum Pr(c) when using counts to set prior

# options to reduce memory use
parser.add_argument("--nSnpPerRoundBF", type=int, default=500000) # how many snps to get BF for at a time

# output options
parser.add_argument("--verbose", type=int, default=0)


args = parser.parse_args()
if args.verbose:
    print(vars(args),file=sys.stderr)
    print("\n", file=sys.stderr)

if args.useSimulatedData:
    if args.nTraits == None:
        msg = "Must specify number of traits to simulate"
        raise argparse.ArgumentTypeError(msg)

if args.useSimulatedData == False:
    args.nTraits = len(args.traits)
    if args.Sigma_Z_file == None:
        # by default, write sample covariance using exptName
        args.Sigma_Z_file = args.outdir+"/"+args.exptName+".Sigma_Z.txt"
    if args.traitPath == None:
        msg = "Must specify the directory containing GWAS summary statistics with --traitPath"
        raise argparse.ArgumentTypeError(msg)
    if args.gwasFormat == None:
        msg = "Must specify format of GWAS summary statistics"
        raise argparse.ArgumentTypeError(msg)



###### simulate or read gwas data #######

configs = list(itertools.product((0,1), repeat=args.nTraits))
altConfigs = list(configs)
altConfigs.remove((0,)*args.nTraits) # list of tup
altConfigsArr = np.array(altConfigs) # np arr
nAltConfigs = len(altConfigs)

# simulate data
if args.useSimulatedData:
    # if no configuration priors specified, set to 0.005 for each alt
    if args.truePriorFile == None:
        truePriorsD = {} # ex/ truePriorsD[(0,-1,1)]=0.01
        for c in configs:
            truePriorsD[c] = 0.005
        truePriorsD[(0,)*args.nTraits] = 1 - (0.005*(2**args.nTraits-1)) # binary traits
        print("Using default simulation priors p(c) of 0.005",file=sys.stderr)

    else:
        truePriorsD = readWtsFile(args.truePriorFile, args.nTraits)

    # also read Sigma_Z from file (if none specified, will be identity)
    Sigma_Z = readSigma_Z(args.Sigma_Z_file, args.nTraits)

    gwasZscores, gwasPvals, trueConfigs = simulateGwasData(args.nSnp, args.nTraits, args.sigmasq_mu_sim, Sigma_Z, truePriorsD, configs, t1trueSigmasq = args.t1trueSigmasq)
    print("Simulated GWAS z-scores and pvals", file=sys.stderr)
  

else:
    snpIDs, gwasZscores, gwasPvals, Sigma_Z = readInGwasData(args.nSnp,
       args.traits, args.traitPath, args.Sigma_Z_file, args.gwasFormat) # also writes estimated Sigma_Z to specified file
    if args.verbose:
        print("Read GWAS z-scores, pvals, and Sigma_Z...", file=sys.stderr)

# get weights, if using initial priors
if args.initialPriorFile:
    weights_A = readWtsFile(args.initialPriorFile, args.nTraits)

else: # use gwas stats to set P(c)
    if args.verbose:
        print("Getting GWAS counts for initial P(c) , over all snps....", file=sys.stderr)
    weights_A = estimateWtsFromGwas(gwasPvals, configs, thresh = args.countPriorThresh, minPrior=args.minPriorC, avg1active = args.avg1active)

print("Computing CONFIT test statistic...", file=sys.stderr)

# to limit memory use, only do up to 500k (or however many specified) snps at a time
BFs_A = np.empty(args.nSnp)
start = 0
nSnps_block = min(args.nSnpPerRoundBF, args.nSnp)
while start < args.nSnp:
    if args.verbose:
        print("Computing CONFIT statistics, on SNP: %d..." % start, file=sys.stderr)
    end = min(start + nSnps_block, args.nSnp)
    gwasZscores_block = gwasZscores[start:end,]
    BFs_A[start:end] = computeCONFITstats(gwasZscores_block, weights_A, args.sigmasq_mu_guess, end-start, args.nTraits, Sigma_Z)
    start += nSnps_block


################# WRITE TO FILE #################

headerL = ['gwas_z'+str(i) for i in range(args.nTraits)] +\
            ['gwas_p'+str(i) for i in range(args.nTraits)] +\
             ['confit_stat']
gwasZscores = [tup2strc(i, d=" ") for i in gwasZscores]
gwasPvals = [tup2strc(i, d=" ") for i in gwasPvals]
BFs_A = [str(i) for i in BFs_A]

if args.useSimulatedData:
    headerL = ['true_c'+str(i) for i in range(args.nTraits)]+headerL
    trueConfigs = [tup2strc(i, d=" ") for i in trueConfigs] # (1,-1) -> 1 -1
    zippedL = zip(trueConfigs, gwasZscores, gwasPvals, BFs_A)

else:
    zippedL = zip(snpIDs, gwasZscores, gwasPvals, BFs_A)

fbase = args.outdir + "/" + args.exptName 
fname = fbase + "_confit.txt"
with open(fname,'w') as outf:
    outf.write(" ".join(headerL) + "\n")
    for i in range(args.nSnp):
        outf.write(" ".join(zippedL[i]) + "\n")

# write weights to file
priorsL = sorted(list(weights_A.iteritems())) # estimated priors

if args.useSimulatedData:
    truePriorsL = sorted(list(truePriorsD.iteritems()) )
    priorsL = zip(truePriorsL, priorsL)
    priorsL = [(" ".join([str(ci) for ci in c]), str(t), str(w)) for ((c,t),(c2,w)) in priorsL]
else:
    priorsL = [(" ".join([str(ci) for ci in c]), str(w)) for c,w in priorsL] 

fname = fbase + "_wt.txt" #_wts.txt"
with open(fname,'w') as outf:
    for i in range(len(priorsL)):
        outf.write(" ".join(priorsL[i]) + "\n")

