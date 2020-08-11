#!/usr/bin/env python 

'''
_version=0.0.1
[Usage]: ./keggenrich.py kegg_annot_file kodb pathway.ko
'''
import os, sys

currsd = os.path.dirname(os.path.realpath(__file__))

def enrich(annotf,kodb,pathway,outfile=None):

	if not outfile:
		outfile = annotf+".kg"
	os.system("%s/keggCountal.py %s %s %s -o %s" % (currsd,annotf,kodb,pathway,outfile))
	sys.stderr.write("countall done\n")
	os.system("%s/mykeggenrich.R %s" % (currsd,outfile))
	sys.stderr.write("enrich R done\n")

def main(argv):
	
	import argparse

	parser = argparse.ArgumentParser(description="kegg enrich")
	parser.add_argument('annotf',help="genes to be enrich. Typically ko annot file of different expression genes")
	parser.add_argument('kodb',nargs='?',help="ko annot file of all genes")
	parser.add_argument('pathway2ko',help="path-ko file")
	parser.add_argument('-o','--outfile',nargs='?',help="output file",default=sys.stdout)
	args = parser.parse_args(argv[1:])
	enrich(args.annotf,args.kodb,args.pathway2ko,args.outfile)

if __name__ == '__main__':

	import sys

	main(sys.argv)
