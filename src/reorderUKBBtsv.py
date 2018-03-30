# re-order all snps by chromosome and position
# (assumes all traits have same SNPs in each file)
# input: path/to/file.tsv

from __future__ import print_function

import numpy as np
import sys

def rewrite(tsvName, outName, orderL):
	# read in unsorted tsv
	i = 0
	tsvD = {}

	with open(tsvName) as f:
		header = f.readline()
		for line in f:
			# get dict of tsvD[snp]: line
			variant, rsid, rest = line.split(None,2)
			tsvD[rsid] = line

			if (i % 5000000) == 0:
				print("reading to rewrite, line %d" % i, file=sys.stderr)
			i += 1

	# write snps in order of orderL
	with open(outName,"w") as outf:
		outf.write(header)
		for (snp, chro, pos) in orderL:
			outf.write(tsvD[snp])
	return

# list of (SNP, chr, pos) [char, int, int]
orderL = []

i = 0
tsvName = sys.argv[1] # get snps from the input TSV, then reorder

with open(tsvName) as f:
	f.readline() # skip header
	for line in f:
		variant, rsid, rest = line.split(None,2)
		chro, pos, rest = variant.split(":",2) #ex/ "1:123456:G:C"
		orderL.append((rsid, int(chro), int(pos)))
		if (i % 1000000) == 0:
			print("reading first file to get list of snps, line %d" % i, file=sys.stderr)
		i += 1

# sort by chr, then by pos
print("sorting by position...", file =sys.stderr)
orderL = sorted(orderL, key = lambda tup: (tup[1], tup[2]))
print("sorted snps in order by position", file =sys.stderr)


outName = tsvName + ".sorted"
rewrite(tsvName, outName, orderL)
print("rewrote %s tsv in order" % trait, file=sys.stderr)


