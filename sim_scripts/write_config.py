#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, getpass, re, argparse
import dendropy
from dendropy.simulate import treesim

parser = argparse.ArgumentParser(description='Simulate constrained coalescent trees for a species tree')

parser.add_argument("-treefile", "--treefilepath", help="Tree file path", required = True)
parser.add_argument("-errorfile", "--errorfilepath", help="Error file path", required = False)
parser.add_argument("-v", "--variablesites", help="Number of variable sites", type=int, action = "store", required = True)
parser.add_argument('-ref', '--referencename', help="Reference base name", required = True)
parser.add_argument('-path', '--referencepath', help="Reference path", required = True)
parser.add_argument('-o', "--outputdir", action="store", help="output file [required]", required=True)
parser.add_argument("-rate", "--ratemat", help="Rate Matrix file", required = True)
parser.add_argument("-cl", "--coverage_low", help="Lower bound of coverage", type=int, action = "store", required = False)
parser.add_argument("-ch", "--coverage_high", help="Upper bound of coverage", type=int, action = "store", required = False)
parser.add_argument("-cs", "--coverage_step", help="Step size for binning coverage", type=int, action = "store", required = False)
parser.add_argument("-g", "--gamma", help="Gamma Shape", type=int, action = "store", required = True)
parser.add_argument("-r", "--read", help="Read length", type=int, action = "store", required = True)
parser.add_argument("-f", "--frag", help="Fragment Size", type=int, action = "store", required = True)
parser.add_argument("-s", "--stddevfrag", help="Stdev Fragment Size", type=int, action = "store", required = True)
parser.add_argument("-p", "--platform", help="Read platform (options: MSv3, GA1, GA2, HS10, HS20, HS25, HSXn, HSXt, MinS, MSv1, MSv2, MSv3, NS50)", action = "store", required = True)
parser.add_argument("-pre", "--prefix", help="Prefix of Simulated Data", action = "store", required = True)

args = parser.parse_args()

# Set the default values:
tree = args.treefilepath
error = args.errorfilepath
variable = args.variablesites
reference = args.referencename
refpath = args.referencepath
output = args.outputdir
coverage_low = args.coverage_low
coverage_high = args.coverage_high
coverage_step = args.coverage_step
gamma = args.gamma
readlength = args.read
fragment = args.frag
stddev = args.stddevfrag
platform = args.platform
prefix = args.prefix

print "treefile_path =", tree
print "number_of_variable_sites =", variable
print "base_genome_name =", reference
print "base_genome_path =", refpath
print "output_dir =" , output
with open(args.ratemat, 'r') as rate:
        print "rate_matrix =", rate.read()
print "coverage_low =", coverage_low
print "coverage_high =", coverage_high
print "coverage_step =", coverage_step
print "prefix =", prefix
#print "gamma_shape =", gamma
print "error_model1 =", error
print "error_model2 =", error
print "read_length =", readlength
print "fragment_size =", fragment
print "stdev_frag_size =", stddev
print "readProfile =", platform
