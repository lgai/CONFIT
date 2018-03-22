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
parser = argparse.ArgumentParser(description="powersim")

parser.add_argument("--exptName", type=str, required=True)

# simulation options
parser.add_argument("--nTraits",  type=int, default=2)

# TODO fix option so that they're True/False
parser.add_argument("--useSimulatedData", type=str, default='False')

parser.add_argument("--nSnp",  type=int, default=5*(10**5))
parser.add_argument("--mu_sim_prior", type=str, default="_CE") # "_USE_CONSTANT" deprecated
parser.add_argument("--mu_sim_const", type=float, default=5.2) # deprecated
parser.add_argument("--sigmasq_mu_sim", type=float, default=25) # should be 25 by default
parser.add_argument("--t1trueSigmasq", type=float, default=None) # if t1 should have diff simulation sigmasq than other traits, set it here

parser.add_argument("--Sigma_Z_file", type=str, default=None) # name of Sigma_ to use for simulation, or where to write Sigma_Z estimated from Finland data. If none and simulating data, use identity
parser.add_argument("--truePriorFile", default=None)


parser.add_argument("--traits", type=str, nargs='*', default=None) # user must edit path name
parser.add_argument("--gwasFormat", type=str, default=None) # default is pylmm, can also use UKbiobank_neale

# options for CONFIT test statistic
parser.add_argument("--mu_guess_prior", type=str, default="_USE_CONSTANT")
parser.add_argument("--mu_guess_const", type=float, default=5.2) # deprecated
parser.add_argument("--sigmasq_mu_guess", type=float, default=None)

# options for prior Pr(c) estimation
parser.add_argument("--initialPriorFile", type=str, default=None)
parser.add_argument("--pAlt_init", type=float, default=0.001) # 1 - Pr_0(c_0) # deprecated
parser.add_argument("--numEstPriorIterations",type=int,default=0) # how many times to iterate pr(c)
parser.add_argument("--useSigOnlyIterative",type=int,default=0) # whether to use only significant snps when doing iteration
parser.add_argument("--useCountPrior",type=int,default=0) # whether to use gwas count prior. default is to split by number of active traits
parser.add_argument("--countPriorThresh",type=float,default=1e-4) # threshold for count prior
parser.add_argument("--priorAll1sOnly", type=int, default=0) # whether to only use prior where the only alternate model is all 1's

parser.add_argument("--avg1active",type=int,default=0) # whether to average over configs w 1 active 
parser.add_argument("--avgNactive",type=int,default=0) # whether to average over singles, pairs, etc

parser.add_argument("--gwasWtThresh", type=float, default=10**(-4))  # gwas pval thresh for estimating priors Pr(c)
parser.add_argument("--minPrior", type=float, default=1e-8) # deprecated

# options to reduce memory use
parser.add_argument("--nSnpPerRoundBF", type=int, default=500000) # how many snps to get BF for at a time

# output options
parser.add_argument("--outdir",type=str,required=True)

args = parser.parse_args()
print(vars(args),file=sys.stderr)
print("\n", file=sys.stderr)

# prob could use as flag instead
if args.useSimulatedData == 'True':
    args.useSimulatedData == True
else:
    args.useSimulatedData = False
    args.nTraits = len(args.traits)
    args.nSnp = 331476 # hardcoded number of snps in real data instead of counting each time
    if args.gwasFormat == "UKbiobank_neale":
        args.nSnp = 10894596 

if args.useSimulatedData == False and args.Sigma_Z_file == None:
    msg = "Must specify Sigma_Z_file location if using real data"
    raise argparse.ArgumentTypeError(msg)

######## READ/GET GWAS DATA, COMPUTE PRIORS AND TEST STAT #########

# binary traits
configs = list(itertools.product((0,1), repeat=args.nTraits))
altConfigs = list(configs)
altConfigs.remove((0,)*args.nTraits) # list of tup
altConfigsArr = np.array(altConfigs) # np arr
nAltConfigs = len(altConfigs)

# simulate data
if args.useSimulatedData:
    # if no config priors specified, set to 0.005 for each alt
    if args.truePriorFile == None:
        truePriorsD = {} # ex/ truePriorsD[(0,-1,1)]=0.01
        for c in configs:
            truePriorsD[c] = 0.005
        truePriorsD[(0,)*args.nTraits] = 1 - (0.005*(2**args.nTraits-1)) # binary traits
        print("using default simulation priors p(c) of 0.005",file=sys.stderr)

    else:
        truePriorsD = readWtsFile(args.truePriorFile, args.nTraits)

    # for simulated data, Sigma_Z is identity by default
    Sigma_Z = readSigma_Z(args.Sigma_Z_file, args.nTraits)

    gwasZscores, gwasPvals, trueConfigs = simulateGwasData(args.nSnp, args.nTraits, args.mu_sim_prior, args.mu_sim_const, args.sigmasq_mu_sim, Sigma_Z, truePriorsD, configs, t1trueSigmasq = args.t1trueSigmasq)
    print("Simulated GWAS z-scores and pvals")
  

else:
    print("Reading GWAS data...")
    gwasZscores, gwasPvals, Sigma_Z = readInGwasData(args.nSnp, args.traits, args.Sigma_Z_file, gwasFormat = args.gwasFormat, outdir=args.outdir, exptName=args.exptName) 
    # also writes estimated Sigma_Z to file
    print("Got GWAS z-scores, pvals, and Sigma_Z...")
    
# get weights, using initial priors
if args.initialPriorFile:
    initialWtsD = readWtsFile(args.initialPriorFile, args.nTraits)

elif args.priorAll1sOnly == 1:
    # first get counts
    initialWtsD = estimateWtsGwas_all(gwasPvals, configs, thresh = args.countPriorThresh, minPrior=1e-8, avg1active = 0, avgNactive = 0) 
    # then set P(all 1's ) = 1 - P(c0), rest to 0
    for c in altConfigs:
        initialWtsD[c] = 0.0
    initialWtsD[(1,)*args.nTraits] = 1 - initialWtsD[(0,)*args.nTraits]
    print("Setting P(all 1) = %f, rest = 0" % initialWtsD[(1,)*args.nTraits], file=sys.stderr)

elif args.useCountPrior: # use counts as initial wts (or final wts if num iter = 0)
    print("Getting GWAS counts for initial P(c) , over all snps....", file=sys.stderr)
    initialWtsD = estimateWtsGwas_all(gwasPvals, configs, thresh = args.countPriorThresh, minPrior=1e-8, avg1active = args.avg1active, avgNactive = args.avgNactive)  # minPrior is deprecated

else:
    initialWtsD = makeInitialPriorsD(args.nTraits, pAlt_init = args.pAlt_init)
    print("Using makeInitialPriorsD() to make default prior....", file=sys.stderr)

weights_A = initialWtsD

# iteratively update priors using entire genome 
for i in range(args.numEstPriorIterations):
    print("estimating wts, iter %d" % i, file=sys.stderr)
    weights_A = estimateWtsSoft_CE(gwasZscores, gwasPvals, weights_A, args.sigmasq_mu_guess, args.nSnp, args.nTraits, Sigma_Z, useSigOnly = args.useSigOnlyIterative)

# get BFs
if args.mu_guess_prior == "_USE_CONSTANT": #  deprecated
    BFs_A = computeBFs_mu_guess_const(gwasZscores, weights_A, args.mu_guess_const, args.nSnp, args.nTraits, Sigma_Z) # use specified Sigma_Z instead of estimating
else: # assume prior for mu_guess
    # to limit memory use, only do up to 500k snps at a time
    BFs_A = np.empty(args.nSnp)
    start = 0
    nSnps_block = min(args.nSnpPerRoundBF, args.nSnp)
    while start < args.nSnp:
        print("computing BFs, start of block: %d..." % start, file=sys.stderr)
        end = min(start + nSnps_block, args.nSnp)
        gwasZscores_block = gwasZscores[start:end,]
        BFs_A[start:end] = computeBFs_mu_guess_prior(gwasZscores_block, weights_A, args.sigmasq_mu_guess, end-start, args.nTraits, Sigma_Z)
        start += nSnps_block


################# WRITE TO FILE #################

headerL = ['gwas_z'+str(i) for i in range(args.nTraits)] +\
            ['gwas_p'+str(i) for i in range(args.nTraits)] +\
             ['BF']
gwasZscores = [tup2strc(i, d=" ") for i in gwasZscores]
gwasPvals = [tup2strc(i, d=" ") for i in gwasPvals]
BFs_A = [str(i) for i in BFs_A]

if args.useSimulatedData:
    headerL = ['true_c'+str(i) for i in range(args.nTraits)]+headerL
    trueConfigs = [tup2strc(i, d=" ") for i in trueConfigs] # (1,-1) -> 1 -1
    zippedL = zip(trueConfigs, gwasZscores, gwasPvals, BFs_A)

else:
    zippedL = zip(gwasZscores, gwasPvals, BFs_A)

fbase = args.outdir + "/" + args.exptName 
fname = fbase + "_BFs.txt"
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

fname = fbase + "_wtsoft.txt" #_wts.txt"
with open(fname,'w') as outf:
    for i in range(len(priorsL)):
        outf.write(" ".join(priorsL[i]) + "\n")

