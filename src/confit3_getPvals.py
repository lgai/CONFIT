# get top statistics from null simulations. In next step, will compare test statistic against these top statistics to get high resolution p-value.

from __future__ import print_function

import os # check if file exists
import sys
import getopt
import argparse
import numpy as np
import numpy.random
import scipy.stats


parser = argparse.ArgumentParser(description="BFpvals")

parser.add_argument("--exptName", type=str, required=True)
parser.add_argument("--nNullRounds", type=int, default=100) # how many nullsim files to read
parser.add_argument("--nNullTotal", type=int, default=5*10**9) # desired p-value resolution
parser.add_argument("--nulldir", type=str, required=True)
parser.add_argument("--outdir",type=str,required=True)




args = parser.parse_args()
print(vars(args),file=sys.stderr)


nTopNull = 1000 
topNullsSoFar = [0]*nTopNull

# maintain that list so it's always sorted
# then take list[0:250] so it's always 250 (or however) long

nNullSoFar = 0
for i in range(1, args.nNullRounds+1): # taskID starts from 1
    nullFile = "%s/nullsim_%s_%d.txt" % (args.nulldir, args.exptName, i)

    if not os.path.isfile(nullFile): # quit if there are null files that failed
        print("Exiting on round %d, file '%s' not found..." % (i, nullFile), file=sys.stderr)
        exit(1) 

    with open(nullFile) as nullf:
        for line in nullf:
            BFnull = float(line.rstrip())
            if BFnull > topNullsSoFar[nTopNull-1]: # if larger than smallest of top nulls so far
                topNullsSoFar.append(BFnull)

            nNullSoFar += 1
            if nNullSoFar == args.nNullTotal:
                break
    
    topNullsSoFar.sort(reverse=True) # sort list after each file, top BF first
    topNullsSoFar = topNullsSoFar[0:nTopNull] # take top 250 or w/e each time
    if i%20 == 0:
        print("reading null CONFIT statistics, round %d" % i, file=sys.stderr)
    if nNullSoFar == args.nNullTotal:
        break 

topNullsSoFar.sort() # now sort so smallest is first

print("read %d null SNPs total" % nNullSoFar, file=sys.stderr)
outName = "%s/%s_topFsnull.txt" % (args.outdir, args.exptName)
with open(outName, "w") as outf:
    outf.write("Top%dNullFs_outOf_%d\n" % (nTopNull,nNullSoFar))

    for BFnull in topNullsSoFar:
        outf.write("%f\n" % BFnull)







