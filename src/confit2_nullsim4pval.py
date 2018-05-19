# generate SNPs under the null

from __future__ import print_function

import sys
import getopt
import argparse
import itertools
import collections
import numpy as np
import numpy.random
import scipy.stats

from confitFunctions import *


# needs python/2.7

#################################################
parser = argparse.ArgumentParser(description="nullsim")
parser.add_argument("--exptName", type=str, required=True)
parser.add_argument("--outdir",  type=str, required=True) # where the output of panel is located - used to get default location of Signa_Z_est and wts_A
parser.add_argument("--nulldir",  type=str, required=True) # where to require output of null simulations
parser.add_argument("--taskID", type=str, default="") # label for each output null file

# can either pass multiple trait names, or specify nTraits
parser.add_argument("--nTraits",  type=int, default=2)
parser.add_argument("--traits", nargs='*', default=None) 
parser.add_argument("--Sigma_Z_file", type=str, default=None) # which Sigma_Z file to use for null simulation. If none, use filename based on exptName. If "identity", use identity matrix
parser.add_argument("--wtsFile", type=str, default=None) # the estimated P(c) file.

parser.add_argument("--sigmasq_mu_guess", type=float, default=25)
parser.add_argument("--fracAssumedGenetic", type=float, default=0.5) # what fraction of cov between traits is assumed genetic

parser.add_argument("--nNullPerRound", type=int, default=2*10**6) # how many snp per iteration of the for loop
parser.add_argument("--nNullRoundsPerJob", type=int, default=25) # how many loops per job

parser.add_argument("--verbose", type=int, default = 1)


args = parser.parse_args()

if args.traits == None and args.nTraits == None:
    msg = "Must specify list of traits OR number of simulated traits"
    raise argparse.ArgumentTypeError(msg)

if args.traits != None:
    nTraits = len(args.traits)
else:
    nTraits = args.nTraits

if args.verbose:
    print(vars(args),file=sys.stderr)

if args.wtsFile == None: # by default, use name from exptName
    args.wtsFile = args.outdir+"/"+args.exptName+"_wt.txt"

estWtsFile = args.wtsFile
taskID = args.taskID

# read Sigma_Z from file
if args.Sigma_Z_file == None: # by default, use name from exptName
    args.Sigma_Z_file = args.outdir+"/"+args.exptName+".Sigma_Z_est.txt"
Sigma_Z = readSigma_Z(args.Sigma_Z_file, nTraits)

weights_A = readWtsFile(estWtsFile, nTraits)
nullFile = "%s/nullsim_%s_%s.txt" % (args.nulldir, args.exptName, args.taskID)

with open(nullFile, 'w') as nullf: 
    for r in xrange(args.nNullRoundsPerJob):
        if args.verbose and r % 10 == 0:
            print("nullsim, round %d..." % r, file=sys.stderr)

        # MI GWAS test stat (abs max z-score across traits)
        # and CONFIT test stat
        MIsAndBFs_null = genBF_nullsim(weights_A, args.nNullPerRound, nTraits, Sigma_Z, args.sigmasq_mu_guess, args.fracAssumedGenetic)

        for MInull, BFnull in MIsAndBFs_null:
            nullf.write(str(MInull) + "\t" + str(BFnull) + "\n")

if args.verbose:
    print("done", file=sys.stderr)

