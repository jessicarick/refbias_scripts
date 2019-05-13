###Write python script to create partition file for raxml - Use intervals
#! /usr/bin/env python
# -*- coding: utf-8 -*-
#usage python snp_part_raxml.py -snpfile snpfile -namesfile namesfile

import sys, getpass, re, argparse, numpy, scipy

parser = argparse.ArgumentParser(description='Write Raxml Partition File')

parser.add_argument("-snpfile", "--snpfilepath", help="SNP file path", required = True)
parser.add_argument("-namesfile", "--namesfilepath", help="Locus Names file path", required = True)

args = parser.parse_args()

# Set the default values:
snpfile = args.snpfilepath
names = args.namesfilepath

with open(snpfile, 'r') as f:
    polyShape = []
    for line in f:
        line = line.split() # to deal with blank 
        if line:            # lines (ie skip them)
            line = [int(i) for i in line]
            polyShape.append(line)

with open(names, 'r') as f: locnames = f.read().split()
            
snps = numpy.array(line)
loci = numpy.cumsum(line)

for i in range(0,snps.size): 

	if i == 0:
		start = 1 
		end = loci[i] 
		print "DNA,", locnames[i], ; print "=", start, ; sys.stdout.softspace=False; print "-", ; sys.stdout.softspace=False; print end
	elif i != 0:
		start = loci[i-1] + 1 
		end = loci[i] 
		print "DNA,", locnames[i], ; print "=", start, ; sys.stdout.softspace=False; print "-", ; sys.stdout.softspace=False; print end
