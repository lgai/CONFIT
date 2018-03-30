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
parser.add_argument("--nulldir",type=str,required=True)
parser.add_argument("--nNullLowRes", type=int, default=10**8
) # the max res for low resolution p-values

args = parser.parse_args()
print(vars(args),file=sys.stderr)
print("\n",file=sys.stderr)
nNullLowRes = args.nNullLowRes

if args.outdir[-1] = "/":
    args.outdir = args.outdir[:-1] # trim trailing slash

# read the BFs from file to get BFs_A
BFs_A = []
with open(args.outdir + "/" + args.exptName + "_confit.txt") as f:
    f.readline() # skip header
    for line in f:
        BFs_A.append(float(line.rstrip().split()[-1]))

# read the top Null BFs from file in ascending order 
topBFsNull = []
with open("%s/%s_topFsnull.txt" % (args.outdir, args.exptName)) as f:
    f.readline() # skip header
    for line in f:
        topBFsNull.append(float(line.rstrip()) )
nTopNull = len(topBFsNull)


BFs_A = np.array(BFs_A)
nSnp = len(BFs_A)
BFs_A = BFs_A.reshape( (1,nSnp) )


## compute p-value with 2e-10 resolution using top null values
# 500k at once to limit memory use
start = 0
nSnps_block = min(500000, nSnp)
nGreaterTop = np.empty(nSnp)
print("start %d, nSnp %d" % (start,nSnp) )
while start < nSnp:
    print("Getting high res pval for snp %d ..." % start, file=sys.stderr)
    end = min(start + nSnps_block, nSnp)
    BFs_A_block = BFs_A[0,start:end]
    topBFsNull = np.array(topBFsNull).reshape((nTopNull,1))
    nGreaterTop[start:end] = (np.sum(BFs_A_block <= topBFsNull, axis = 0)).reshape(end-start) # shape is (1, nSnp_block)
    start += nSnps_block


print("Got pvals against %d top null values " % nTopNull, file=sys.stderr)


### compute p-values with less precision for less significant SNPs 

# read null BFs from first few rounds of null sim to get low res pval
nullBFs = np.empty((nNullLowRes))
i = 0

for jobNo in range(1,500): # look through up to 500 null files, until you have 1e8 or w/e pvals
    if i == nNullLowRes:
        break
    with open("%s/nullsim_%s_%d.txt" % (args.nulldir, args.exptName, jobNo)) as f:
        for line in f:
            if i == nNullLowRes:
                break
            nullBFs[i] = float(line.rstrip())
            i += 1
            if i % 5*10**7 == 0:
                print("Read %d snps for low res so far..." % i,file=sys.stderr)



print("Read %d snps for low res, stopping" % i,file=sys.stderr)

nGreaterLowRes = np.ones((nSnp))

for i in range(nSnp):
    if i%500000==0:
        print("Getting low res pval for snp %d..." % i, file=sys.stderr)

    # compare each snp to subset of null values

    BF_i = BFs_A[0,i]

    nGreaterLowRes = np.sum(BF_i <= nullBFs[0:10**3])
    if nGreaterLowRes_X > 10: #  if pval > 10/1000, use 1e3 resolution
        nGreaterLowRes[i] = nGreaterLowRes_X * int(nNullLowRes/10**3)
        continue
    
    nGreaterLowRes_X = np.sum(BF_i <= nullBFs[0:10**4])
    if nGreaterLowRes_X > 10:
        nGreaterLowRes[i] = nGreaterLowRes_X * int(nNullLowRes/10**4)
        continue

    nGreaterLowRes_X = np.sum(BF_i <= nullBFs[0:10**5])
    if nGreaterLowRes_X > 20:
        nGreaterLowRes[i] = nGreaterLowRes_X * int(nNullLowRes/10**5)
        continue

    nGreaterLowRes_X = np.sum(BF_i <= nullBFs[0:10**7])

    if nGreaterLowRes_X > 100: 
        nGreaterLowRes[i] = nGreaterLowRes_Xi * int(nNullLowRes/10**7)
        continue
    
    nGreaterLowRes[i] = np.sum(BF_i <= nullBFs) 



# combine with the _BFs.txt file
with open(args.outdir + "/" + args.exptName + "_confitwpvals.txt", "w") as outf:
    with open(args.outdir + "/" + args.exptName + "_confit.txt") as f:
        header = f.readline().rstrip()
        header += " nNullGreaterThan_top nNullGreaterThan_lowres nNullGreaterThan_lowres_flat\n"
        outf.write(header)
        i = 0 # to get corresponding pval
        for line in f: 
            outf.write("%s %d %d %d\n" % (line.rstrip(), nGreaterTop[i], nGreaterLowRes[i], nGreaterLowRes_flat[i]))
            i += 1

