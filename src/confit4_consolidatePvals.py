from __future__ import print_function

import os # to check if file exists
import sys
import getopt
import argparse
import numpy as np
import subprocess # to execute shell command

parser = argparse.ArgumentParser(description="BFpvals_consolidate")

parser.add_argument("--exptName", type=str, required=True)
# need <= 5 out of 10^8 to be gwas sig, or 2 out of 4*10^7
parser.add_argument("--outdir",type=str,required=True)
parser.add_argument("--nulldir",type=str,required=True) # where the null files are located
parser.add_argument("--useSimulatedData",type=int,default=0) # whether or not simulated data was used (affects file format)
parser.add_argument("--nTraits",type=int,required=True) # how many traits (affects file format)
parser.add_argument("--nNullLowRes", type=int, default=10**8
) # the max res for low resolution p-values
parser.add_argument("--nNullFullRes", type=int, default=5*10**9) # for high resolution p-values

args = parser.parse_args()
print(vars(args),file=sys.stderr)
print("\n",file=sys.stderr)
nNullLowRes = args.nNullLowRes
nNullFullRes = args.nNullFullRes

nTraits = args.nTraits

if args.outdir[-1] == "/":
    args.outdir = args.outdir[:-1] # trim trailing slash

# read the gwas max abs z-score and BFs from file
# of the test snps
MIs_A = []
BFs_A = []
with open(args.outdir + "/" + args.exptName + "_confit.txt") as f:
    f.readline() # skip header
    for line in f:
        lineL = line.rstrip().split()
        # get max abs gwas z-score across traits
        if args.useSimulatedData: # first nTraits columns are true config
            absMI_zs = [abs(float(z)) for z in lineL[nTraits:nTraits*3:2]]
        else: # first col is rsID
            absMI_zs = [abs(float(z)) for z in lineL[1:nTraits*2+1:2]]
        print(absMI_zs) # TODO
        MIs_A.append(max(absMI_zs))
        BFs_A.append(float(lineL[-1]))


# read top null MI GWAS and BFs from file in ascending order
topMIsNull = []
topBFsNull = []
with open("%s/%s_topMIGWASnull.txt" % (args.outdir, args.exptName)) as f:
    f.readline() # skip header
    for line in f:
        topMIsNull.append(float(line.rstrip()) )
nTopNull = len(topMIsNull)
topMIsNull = np.array(topMIsNull).reshape((nTopNull,1))

topBFsNull = []
with open("%s/%s_topCONFITnull.txt" % (args.outdir, args.exptName)) as f:
    f.readline() # skip header
    for line in f:
        topBFsNull.append(float(line.rstrip()) )
nTopNull = len(topBFsNull)
topBFsNull = np.array(topBFsNull).reshape((nTopNull,1))

nSnp = len(MIs_A)
MIs_A = np.array(MIs_A)
MIs_A = MIs_A.reshape ( (1, nSnp) )
BFs_A = np.array(BFs_A)
BFs_A = BFs_A.reshape( (1,nSnp) )

## compute high resolution p-value using top null values
# 500k at once to limit memory use (can change)
start = 0
nSnps_block = min(500000, nSnp)
nGreaterTopMI = np.empty(nSnp)
nGreaterTopBF = np.empty(nSnp)

print("start %d, nSnp %d" % (start,nSnp) )
while start < nSnp:
    print("Getting high res pval for snp %d ..." % start, file=sys.stderr)
    end = min(start + nSnps_block, nSnp)

    MIs_A_block = MIs_A[0,start:end]
    nGreaterTopMI[start:end] = (np.sum(MIs_A_block <= topMIsNull, axis = 0)).reshape(end-start) # shape is (1, nSnp_block)

    BFs_A_block = BFs_A[0,start:end]

    nGreaterTopBF[start:end] = (np.sum(BFs_A_block <= topBFsNull, axis = 0)).reshape(end-start) # shape is (1, nSnp_block)
    
    start += nSnps_block

pvalsHighResMI = 1.0/nNullFullRes*nGreaterTopMI
pvalsHighResBF = 1.0/nNullFullRes*nGreaterTopBF

print("Got pvals against %d top null values " % nTopNull, file=sys.stderr)


### compute p-values with less precision for less significant SNPs 

# read null test stats from first few rounds of null sim to get low res pval
nullMIs = np.empty((nNullLowRes))
nullBFs = np.empty((nNullLowRes))
i = 0

for jobNo in range(1,500): # look through up to 1000 null files (can change), until you have enough null values for low resolution
    if i == nNullLowRes:
        break
    with open("%s/nullsim_%s_%d.txt" % (args.nulldir, args.exptName, jobNo)) as f:
        for line in f:
            if i == nNullLowRes:
                break
            lineL = line.rstrip().split()
            nullMIs[i] = float(lineL[0])
            nullBFs[i] = float(lineL[-1])
            i += 1
            if i % 5*10**7 == 0:
                print("Read %d snps for low res so far..." % i,file=sys.stderr)

print("Read %d snps for low res, stopping" % i,file=sys.stderr)

 
def getLowResPval(testStat, nullStats, nNullLowRes):
    # testStat should be a float
    # nullStats should be (nNullLowRes,)

    nSampleL = [10**3, 10**4, 10**5, 10**6, 10**7]

    for nSample in nSampleL:
        if nSample >= nNullLowRes:
            break
        nGreater_sm = np.sum(testStat <= nullStats[0:nSample])
        if nGreater_sm > 0.005*nSample: # e.g. if pval > 10/1000
            nGreater = nGreater_sm*int(nNullLowRes/nSample)
            return 1.0*nGreater/nNullLowRes

    # if pval is small/haven't returned yet, use max nNullLowRes
    nGreater = np.sum(testStat <= nullStats[0:nNullLowRes])

    return 1.0*nGreater/nNullLowRes


# compare each snp to a sample of null values to get low res pvals

pvalsLowResMI = np.ones((nSnp))
pvalsLowResBF = np.ones((nSnp))

for i in range(nSnp):
    if i%500000==0:
        print("Getting low res pval for snp %d..." % i, file=sys.stderr)
   
    MI_i = MIs_A[0,i]
    pvalsLowResMI[i] = getLowResPval(MI_i, nullMIs, args.nNullLowRes)

    BF_i = BFs_A[0,i]
    pvalsLowResBF[i] = getLowResPval(BF_i, nullBFs, args.nNullLowRes)

# write statistics with pvals (for both MI GWAS and CONFIT)
with open(args.outdir + "/" + args.exptName + "_confitwpvals.txt", "w") as outf:
    with open(args.outdir + "/" + args.exptName + "_confit.txt") as f:
        header = f.readline().rstrip()
        header += " pval_MIGWAS_lowres pval_MIGWAS_highres pval_CONFIT_lowres pval_CONFIT_highres\n"
        outf.write(header)
        i = 0 # to get corresponding pval
        for line in f: 
            outf.write("%s %s %s %s %s\n" % (
                line.rstrip(), 
                str(pvalsLowResMI[i]),
                str(pvalsHighResMI[i]),
                str(pvalsLowResBF[i]),
                str(pvalsHighResBF[i])) )
            i += 1

print("done", file=sys.stderr)

