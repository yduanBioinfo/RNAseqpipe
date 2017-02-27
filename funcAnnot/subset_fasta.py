#!/usr/bin/env python

import sys

'''infiles fasta_file IDsfile
   outfile fasta_file contains all IDs and its sequence in IDsfile
   usage
   ./subset_fasta.py fasta_file IDsfile outfile
'''

infasta = open(sys.argv[1])
IDfile = open(sys.argv[2])
outfasta = open(sys.argv[3],'w')

myIDs = [i.strip() for i in IDfile]
writefile = False

for eachline in infasta:
	if eachline[0] == ">":
		if eachline.strip().lstrip(">") in myIDs:
			writefile = True
		else:
			writefile = False
	if writefile:
		outfasta.write(eachline)

infasta.close()
outfasta.close()