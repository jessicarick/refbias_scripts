#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys, getpass, re, argparse
import dendropy
from dendropy.simulate import treesim

parser = argparse.ArgumentParser(description='Simulate contained coalescent trees for a species tree')

parser.add_argument("-treefile", "--treefilepath", help="Tree file path", required = True)
parser.add_argument("-errorfile", "--errorfilepath", help="Error file path", required = False)
parser.add_argument("-v", "--variablesites", help="Number of variable sites", type=int, action = "store", required = True)
parser.add_argument('-ref', '--referencename', help="Reference base name", required = True)
parser.add_argument('-path', '--referencepath', help="Reference path", required = True)
parser.add_argument('-o', "--outputdir", action="store", help="output file [required]", required=True)
parser.add_argument("-rate", "--ratemat", help="Rate Matrix file", required = True)
parser.add_argument("-c", "--coverage", help="Coverage", type=int, action = "store", required = True)
parser.add_argument("-g", "--gamma", help="Gamma Shape", type=int, action = "store", required = True)
parser.add_argument("-r", "--read", help="Read length", type=int, action = "store", required = True)
parser.add_argument("-f", "--frag", help="Fragment Size", type=int, action = "store", required = True)
parser.add_argument("-s", "--stddevfrag", help="Stdev Fragment Size", type=int, action = "store", required = True)
parser.add_argument("-pre", "--prefix", help="Prefix of Simulated Data", action = "store", required = True)

args = parser.parse_args()

# Set the default values:
tree = args.treefilepath
error = args.errorfilepath
variable = args.variablesites
reference = args.referencename
refpath = args.referencepath
output = args.outputdir
coverage = args.coverage
gamma = args.gamma
readlength = args.read
fragment = args.frag
stddev = args.stddevfrag
prefix = args.prefix

print "treefile_path =", tree
print "number_of_variable_sites =", variable
print "base_genome_name =", reference
print "base_genome_path =", refpath
print "output_dir =" , output
with open(args.ratemat, 'r') as rate:
        print "rate_matrix =", rate.read()
print "coverage =", coverage
print "prefix =", prefix
print "gamma_shape =", gamma
print "error_model1 =", error
print "error_model2 =", error
print "read_length =", readlength
print "fragment_size =", fragment
print "stdev_frag_size =", stddev
