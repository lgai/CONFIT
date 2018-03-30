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
parser.add_argument("--taskID", type=str, default="")

# can either pass multiple trait names, or specify nTraits
parser.add_argument("--nTraits",  type=int, default=2)
parser.add_argument("--traits", nargs='*', default=None) 
parser.add_argument("--Sigma_Z_file", type=str, default=None) # which finland Sigma_Z to use for null simulation. If none, use identity

parser.add_argument("--sigmasq_mu_guess", type=float, default=None)

parser.add_argument("--nNullPerRound", type=int, default=2*10**6) # how many snp per iteration of the for loop
parser.add_argument("--nNullRoundsPerJob", type=int, default=25) # how many loops per job
parser.add_argument("--nulldir",  type=str, required=True) # where to require output of null simulations
parser.add_argument("--wtsFile", type=str, required=True) # the estimated P(c) file 

args = parser.parse_args()

if args.traits == None and args.nTraits == None:
    msg = "Must specify list of traits OR number of simulated traits"
    raise argparse.ArgumentTypeError(msg)

if args.traits != None:
    nTraits = len(args.traits)
else:
    nTraits = args.nTraits

print(vars(args),file=sys.stderr)
print("nTraits %d" % nTraits, file=sys.stderr)


estWtsFile = args.wtsFile % args.exptName


# make 5*10^8 null snps per job, 
taskID = args.taskID


# read Sigma_Z from file, or identity if not specified
Sigma_Z = readSigma_Z(args.Sigma_Z_file, nTraits)

weights_A = readWtsFile(estWtsFile, nTraits)
 

nullFile = "%s/nullsim_%s_%s.txt" % (args.nulldir, args.exptName, args.taskID)

with open(nullFile, 'w') as nullf: 
    for r in xrange(args.nNullRoundsPerJob):
        if r % 10 == 0:
            print("round %d..." % r, file=sys.stderr)

        BFs_null = genBF_nullsim(weights_A, args.nNullPerRound, nTraits, Sigma_Z, args.sigmasq_mu_guess)

        for BFnull in BFs_null:
            nullf.write(str(BFnull) + "\n")

            
