#! /usr/bin/python3
# This script is used to convert sfs by city by habitat obtained by
# ANGSD in EPOS format, for analyse of Ne

import sys

for line in open(sys.argv[1]):
    liste = line.split()
          
out = open(sys.argv[2], "w")
N = -1
for i in liste :
    N=N+1
    if N<=((len(liste)-1)/2) :
            out.write(str(N)+"\t"+str(i.rstrip())+"\n")