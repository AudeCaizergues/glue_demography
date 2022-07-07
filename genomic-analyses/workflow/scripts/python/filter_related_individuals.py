#! /usr/bin/python3
# This script is used to find related individuals after running ngsRelate
# The threshold used here is the r(ab)=0.5 (pairwise relatedness Hedrick et al)
# which represent first degree of relatedness (parent-offspring, or full-sib).
# Run on ngsRelated file output.

import sys
out=open(sys.argv[1]+".related", "w")

for line in open(sys.argv[1],"r") :
	if line.startswith("a	b"):
		out.write(line.rstrip()+"\n")
	if not line.startswith("a	b"):
		values = line.split("\t")
		rab = float(values[14])
		if rab >= 0.5 :
			out.write(line.rstrip()+"\n")
		
