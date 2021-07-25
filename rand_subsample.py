#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, getpass, re, argparse, random, numpy, scipy

parser = argparse.ArgumentParser(description='Subsample n, random, non-overlapping intervals of length l from sequence')

##Define arguments
parser.add_argument("-s", "--start", help="Start of sequence or Range to sample from", type=int, action = "store", required = True)
parser.add_argument("-e", "--end", help="End of sequence or Range to sample from", type=int, action = "store", required = True)
parser.add_argument("-n", "--number", help="Number of samples of length l", type=int, action = "store", required = True)
parser.add_argument('-l', '--length', help="Length of each subsample", type=int, action = "store", required = True)
#parser.add_argument('-o', '--output', dest='o', help="output file [required]", required=True)

args = parser.parse_args()

# Set the default values:
start = args.start
end = args.end
n = args.number
l = args.length
#output = args.output

###Random sampling function using a shift factor to avoid overlap in intervals sampled
def random_subsample(seq, n, l):
    indices = xrange(len(seq) - (l - 1) * n)
    res = []
    shift = 0
    for i in sorted(random.sample(indices, n)):
        i += shift
        res.append(seq[i:i+l])
        shift += l - 1
    return res

   
#print "sequence ='", start, end
#print "Number of samples =", n
#print "Length of Samples =", l

seq = numpy.arange(start, end, 1)

res = random_subsample(seq, n, l)

#print res

for i in range(0, len(res)):

	print str(res[i][0]) + "-" + str(res[i][l-1])

