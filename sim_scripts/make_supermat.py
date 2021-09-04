#!/usr/bin/python

from Bio import AlignIO
import glob
import sys

output = sys.argv[1]

filelist = glob.glob('./gene*.noInv.phy')
i = 0

print("Concatenating gene phylips")

for file in filelist:
#	print(file)
	alignment = AlignIO.read(file, "phylip-relaxed")
	alignment.sort()
#	print("Alignment of length %i" % alignment.get_alignment_length())
	if i==0:
		cat_algn = alignment
	else:
		cat_algn += alignment
	i += 1

print("Concatenated alignment of length %i" % cat_algn.get_alignment_length())

outfh = open(output,"w")
AlignIO.write(cat_algn, outfh, "phylip-relaxed")

outfh.close()
