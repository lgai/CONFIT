from __future__ import print_function

import sys
import getopt
import argparse
import itertools
import collections
import numpy as np
import numpy.random
import scipy.stats
from scipy.stats import multivariate_normal as mvnorm # uses python2.7
import random 


# estimates P(c) 
# by checking where gwas p-val < thresh for each trait, for each snp
def estimateWtsFromGwas(gwasPvals, configs, thresh=1e-4, minPrior=1e-8, avg1active=True):
  
    # get gwas config, convert rows to tuples
    gwasConfigs = -9*np.ones_like(gwasPvals) 

    inds = np.where(gwasPvals > thresh)
    gwasConfigs[inds] = 0
    inds = np.where(gwasPvals <= thresh)
    gwasConfigs[inds] = 1

    gwasConfigs = map(tuple, gwasConfigs)
    
    # count each empirical configuration to get priors
    estWtsD = {}

    countsD = collections.Counter(gwasConfigs)
    nSnpUsed = sum([countsD[c] for c in configs])

    print("Got GWAS counts from all snps, %d snps used" % nSnpUsed, file=sys.stderr)
    for c in configs: # in py3.x, use .items()
        #print("config,count",c,countsD[c],file=sys.stderr)
        estWtsD[c] = max(1.0*countsD[c]/nSnpUsed, minPrior)
    
    # if averaging over configs with single active
    if avg1active:
        print("Averaging weights for config with 1 active trait", file=sys.stderr)
        sum1active = 0
        for c in configs:
            if sum(c)==1:
                sum1active += estWtsD[c]
        for c in configs:
            if sum(c) == 1:
                estWtsD[c] = 1.0*sum1active/len(c)

    # rescale estWtsD so weights sum to 1
    # also set Pr(c0) as frac of snps that are GWAS significant in 1+ trait (or to small value if none are gwas sig)
    altConfigs = list(configs)
    nTraits = len(configs[0])
    altConfigs.remove((0,)*nTraits)
    MIgwasPvals = gwasPvals.min(axis=1)
    snpIndices = np.where(MIgwasPvals < 5e-8)
    nSig = len(snpIndices[0])
    nSnp = len(MIgwasPvals)
    pAlt = max(1.0*nSig/nSnp, minPrior)

    estWtsD[(0,)*nTraits] = 1 - pAlt
    sumGivenSig = sum([estWtsD[c] for c in altConfigs]) # sum p(c|gwas sig) for alt c

    for c in altConfigs:
        estWtsD[c] = 1.0*estWtsD[c]*pAlt/sumGivenSig

    return estWtsD


##############################################################################

# get BFs assuming (z|c) ~ MVN(0, sigmasq*diag(c) + Sigma_e)
# BF = sum_c [P(c) P(z|c,mu_guess)]/P(c0)P(z|c0)
def computeCONFITstats(gwasZscores, estWtsD, sigmasq_mu_guess, nSnp, nTraits, Sigma_Z):
    
    configs = list(itertools.product((0,1),repeat=nTraits)) 
    altConfigs = list(configs)
    altConfigs.remove((0,)*nTraits)
    nBinaryAltConfigs = len(altConfigs)

    # divide off-diagonal elements of empirical covariance by 2 (50% herit)
    Sigma_e = 0.5*Sigma_Z + 0.5*np.identity(nTraits)
    Sigma_e = np.matrix(Sigma_e)

    # get P(z|c) with prior on lambda 
    # assume (z|c_0) ~ N(0, Sigma_e)
    pzcs_null = mvnorm.pdf(gwasZscores, mean=np.zeros(nTraits), cov=Sigma_e)
    pzcs_null = np.array(pzcs_null) # if one at time, must convert float to array

    pzcs_alt = np.zeros((nSnp, nBinaryAltConfigs)) 

    epsilon_sm = 10**-8 # small const for cov matrix

    # for each config, pdf from same mvn for all snp
    # i.e. for each config, pass in 2D (nSNP * nTrait) array of MVN samples to pdf
    # then rejoin pdfs into 2D array of nSnp * nConfig dimension

    for i in range(nBinaryAltConfigs):
        c = altConfigs[i]
        Sigma_c = sigmasq_mu_guess * np.diag(c) + epsilon_sm * np.identity(nTraits)
        Sigma_c = np.matrix(Sigma_c)

        zCov = Sigma_c + Sigma_e # cov for marginal dist
        pzc = mvnorm.pdf(gwasZscores, mean=np.zeros((nTraits,)), cov=zCov)
        pzcs_alt[:,i] = pzc

    
    altPriors = [estWtsD[c] for c in altConfigs] # vector of priors in order of altConfigs
    pzcs_alt_summed = pzcs_alt.dot(altPriors) # get sum of entries p(z,c) = p(c)p(z|c)
    pzcs_null = pzcs_null * estWtsD[(0,)*nTraits] # p(z,c0) = p(c0)p(z|c0)

    for i in range(nSnp):
        if pzcs_null[i] == 0:
            print("At SNP i = %d, p(z,c0) is approx 0, setting p(z,c0) = 1e-100 to avoid divide by 0..." % i, file=sys.stderr)
            pzcs_null[i] = 1e-100 # set as a small positive value to avoid divide by 0

    # Now compute sum_c [P(c) P(z|c,mu_guess)]/P(c0)P(z|c0) 
    BFs = pzcs_alt_summed/pzcs_null

    return BFs

###############################################################

# perform null simulations for significance threshold
# return BFs from 1 panel of null snps in a vector, where z ~ N(0, \Sigma_Z) under null
# later, use null BFs to get CONFIT pval and threshold
def genBF_nullsim(weights_A, nNullPerRound, nTraits, Sigma_Z, sigmasq_mu_guess, returnMI=True):

    gwasZscores = np.empty((nNullPerRound, nTraits))
    for i in range(nNullPerRound):
        gwasZscores[i,:] = np.random.multivariate_normal(np.zeros(nTraits), Sigma_Z) # each z_v is a dim (nTrait) 

    BFs = computeCONFITstats(gwasZscores, weights_A, sigmasq_mu_guess, nNullPerRound, nTraits, Sigma_Z)

    if returnMI: # for MI gwas, return max abs z-score for each snp
        MIs = np.amax(np.abs(gwasZscores),axis=1)
        return(zip(MIs, BFs))
    else:
        return(BFs) 

###########################################################################

## for simulating gwas summary statistics

# get priors from prior file with format
# 0 1 0.01 (e.g. first columns are the configuration, prior is last column
# if prior file also contains simulation priors in 2nd to last col and estimated priors in last column, only use values in last column
def readWtsFile(wtsFile, nTraits):
    weightsD = {}
    with open(wtsFile) as f:
        for line in f.readlines():
            lineL = line.rstrip().split()
            c = tuple([int(ci) for ci in lineL[0:nTraits]])
            prior = float(lineL[-1]) # use estWts if file contains both
            weightsD[c] = prior
    return weightsD

# np doesn't allow sampling from list of tups so have to convert
def tup2strc(config, d=","): # (1,-1) -> "1,-1"
    return d.join([str(ci) for ci in config])

def str2tupc(config, d=","): #"1,-1" -> (1,-1)
    return tuple([int(c) for c in config.split(d)])

def readSigma_Z(Sigma_Z_file, nTraits, verbose=False):
    # read Sigma_Z from file
    if Sigma_Z_file == None: # default is identity
        return np.identity(nTraits)

    Sigma_Z = np.empty((nTraits,nTraits))
    i = 0
    with open(Sigma_Z_file) as f:
        for line in f:
            Sigma_Z[i,:] = [float(Sij) for Sij in line.rstrip().split()]
            i += 1
    if verbose:
        print(Sigma_Z, file=sys.stderr)

    return np.matrix(Sigma_Z)


# using true priors on config, select config for each snp, as tuples
def simulateGwasData(nSnp, nTraits, sigmasq_mu_sim, Sigma_Z, truePriorsD, configs, t1trueSigmasq=None): 
    # using true priors on config, select config for each snp, as tuples
    trueConfigs = numpy.random.choice([tup2strc(c) for c in configs], p=[truePriorsD[c] for c in configs],
                            size=nSnp)
    trueConfigs = map(str2tupc, trueConfigs)

    # draw vectors of NCP from marginal mvn of z
    # (same as drawing mu_sim, then z)
    epsilon_sm = 10**-8
    gwasZscores = np.empty((nSnp, nTraits))
    for i in range(nSnp):
        c = trueConfigs[i]
        Sigma_c = sigmasq_mu_sim * np.diag(c) + epsilon_sm * np.identity(nTraits)
        if t1trueSigmasq and c[0] == 1: # if trait 1 has different variance
            Sigma_c[0,0] = t1trueSigmasq

        zCov = Sigma_c + Sigma_Z # add correlation with Sigma_Z
        gwasZscores[i,:] = np.random.multivariate_normal(np.zeros(nTraits), zCov) # each z_v is a dim (nTrait) vector

    gwasPvals = 2*scipy.stats.norm.cdf(-1.0*np.absolute(gwasZscores))
    return gwasZscores, gwasPvals, trueConfigs



# read summary statistics from formatted file
# for each trait, convert betas and sds to z-scores

def getGwasFileName(trait, traitPath, gwasFormat="pylmm"): # get summary stats file name
    if gwasFormat == "pylmm":
        s = "%s/%s" % (traitPath, trait) 
    elif gwasFormat == "UKBB_sorted": # NOTE: use sorted files so snps in same order
        s = "%s/%s" % (traitPath, trait)
    else:
        print("ERROR: GWAS format '%s' not recognized" % gwasFormat, file=sys.stderr)
        exit(1)
    return s

# get snp ids, z-score (t-score), python pval 
def getSummaryStatsTrait(pheno, traitPath, gwasFormat):
    print("Reading trait %s..." % pheno, file=sys.stderr)
    snpIDs = []
    gwasZscores_t = []
    gwasPvals_t = []
    with open( getGwasFileName(pheno, traitPath, gwasFormat) ) as f:
        f.readline() # skip header
        for line in f:
            lineL = line.rstrip().split()
            if gwasFormat=="pylmm": # default is pylmm
                snpIDs.append(lineL[0])
                gwasZscores_t.append( float(lineL[1])/float(lineL[2]) )
                gwasPvals_t.append( float(lineL[4]) )
            elif gwasFormat=="UKBB_sorted":
                snpIDs.append(lineL[1]) # rs ID is 2nd col
                gwasZscores_t.append( float(lineL[-2]) ) # phesant t-score
                gwasPvals_t.append( float(lineL[-1])) # phesant pval
                
            else:
                print("ERROR: '%s' is not a recognized option for GWAS format" % gwasFormat, file=sys.stderr)
                exit(1)

    #gwasPvals_t = np.array(gwasPvals_t) # can also use p-value computed by input GWAS
    gwasPvals_t = 2*scipy.stats.norm.cdf(-1.0*np.absolute(gwasZscores_t))

    return (snpIDs, gwasZscores_t, gwasPvals_t)


def readInGwasData(nSnp, traitsL, traitPath, Sigma_Z_file, gwasFormat):

    nTraits = len(traitsL)

    # all possible configs
    configs = list(itertools.product((0,1), repeat=nTraits))

    # list of possible configs and alt configs
    altConfigs = list(configs)
    altConfigs.remove((0,)*nTraits) # list of tup
    altConfigsArr = np.array(altConfigs) # np arr
    nAltConfigs = len(altConfigs)

    gwasZscores = np.empty((nSnp, nTraits))
    gwasPvals = np.empty((nSnp, nTraits)) # using pylmm pvals, t-score maybe?

    # get summary stats
    for t in range(nTraits):
        (snpIDs, gwasZscores_t, gwasPvals_t) = getSummaryStatsTrait(traitsL[t], traitPath, gwasFormat)
        gwasZscores[:,t] = gwasZscores_t
        gwasPvals[:,t] = gwasPvals_t

    if Sigma_Z_file == None: # use identity by default
        Sigma_Z_est = np.eye(nTraits)

    # if using summary statistics to estimate Sigma_Z
    else:
        # do corrcoef of nTraits x nSnp array
        gwasZscores_forSigma = gwasZscores
        Sigma_Z_est = np.corrcoef(np.transpose(gwasZscores_forSigma))

    # write estimated Sigma_Z to file - nullsim step will use Sigma_Z
    with open(Sigma_Z_file,"w") as outf:
        for i in range(nTraits):
            outf.write( " ".join([str(Sij) for Sij in  Sigma_Z_est[i,:]]) + "\n")

    Sigma_Z_est = np.matrix(Sigma_Z_est)

    return snpIDs, gwasZscores, gwasPvals, Sigma_Z_est




   