# get top MI GWAS and CONFIT test statistics from null simulations. In next step, will compare test statistic against these top statistics to get high resolution p-value.
# ex/ have top 1000 pvals out of 5e9 null snps. If test stat < (only 3 of the 1000 top values), then test stat < 3 null values out of all 5e9. then pval is P(nullBF < testBF) = 3/5e9, but only had to check against 1000 values instead of 5e9 for each snp

from __future__ import print_function

import os # check if file exists
import sys
import getopt
import argparse
import numpy as np
import numpy.random
import scipy.stats


parser = argparse.ArgumentParser(description="topTestStatistics")

parser.add_argument("--exptName", type=str, required=True)
parser.add_argument("--nulldir", type=str, required=True)
parser.add_argument("--outdir",type=str,required=True)

parser.add_argument("--nNullFiles", type=int, default=100) # how many nullsim files to read
parser.add_argument("--nNullTotal", type=int, default=5*10**9) # desired p-value resolution

parser.add_argument("--nTopNull",type=int,default=1000) # how many of top null



args = parser.parse_args()
print(vars(args),file=sys.stderr)


nTopNull = args.nTopNull
topMINullsSoFar = [0]*nTopNull
topBFNullsSoFar = [0]*nTopNull

# maintain that list so it's always sorted
# then take list[0:250] so it's always 250 (or however) long

nNullSoFar = 0
for i in range(1, args.nNullFiles+1): # taskID starts from 1
    nullFile = "%s/nullsim_%s_%d.txt" % (args.nulldir, args.exptName, i)

    if not os.path.isfile(nullFile): # quit if there are null files that failed
        print("Exiting on file %d, file '%s' not found..." % (i, nullFile), file=sys.stderr)
        exit(1) 

    with open(nullFile) as nullf:
        for line in nullf:
            lineL = line.rstrip().split()
            MInull = float(lineL[0])
            BFnull = float(lineL[1])

            # if larger than smallest top null so far, add it
            if MInull > topBFNullsSoFar[nTopNull-1]: 
                topMINullsSoFar.append(MInull)
            if BFnull > topBFNullsSoFar[nTopNull-1]:
                topBFNullsSoFar.append(BFnull)

            nNullSoFar += 1 # how many null values viewed so far
            if nNullSoFar == args.nNullTotal:
                break
    
    # after each file, take top (largest) abs z-scores or BFs
    topMINullsSoFar.sort(reverse=True) 
    topMINullsSoFar = topMINullsSoFar[0:nTopNull] 
    topBFNullsSoFar.sort(reverse=True) 
    topBFNullsSoFar = topBFNullsSoFar[0:nTopNull] 
    if i % 20 == 0:
        print("reading null CONFIT statistics, round %d" % i, file=sys.stderr)
    if nNullSoFar == args.nNullTotal:
        break 

topMINullsSoFar.sort() # now sort so smallest is first
topBFNullsSoFar.sort()

print("read %d null SNPs total" % nNullSoFar, file=sys.stderr)
outName = "%s/%s_topCONFITnull.txt" % (args.outdir, args.exptName)
with open(outName, "w") as outf:
    outf.write("Top%dNullCONFIT_outOf_%d\n" % (nTopNull,nNullSoFar))
    for BFnull in topBFNullsSoFar:
        outf.write("%s\n" % str(BFnull))

outName = "%s/%s_topMIGWASnull.txt" % (args.outdir, args.exptName)
with open(outName, "w") as outf:
    outf.write("Top%dNullAbsMIGwas_outOf_%d\n" % (nTopNull,nNullSoFar))
    for MInull in topMINullsSoFar:
        outf.write("%s\n" % str(MInull))


print("wrote top MI GWAS and CONFIT test statistics", file=sys.stderr)



