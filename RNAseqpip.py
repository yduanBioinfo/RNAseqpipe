#!/usr/bin/env python
'''
_version = 0.2.3
wrapper for RNAseqpip programs
'''
from __future__ import print_function
import sys, re, os, copy
from threading import Thread
import subprocess

import hisat, cufflinks, NOISeq, htseq
import funcAnnot.b2gprog.GOannot as GOannot
import funcAnnot.kegg.KOannot as KOannot
from funcAnnot.kegg.annot2go_stat import annot2go_stat
from funcAnnot.annot_from_db import annot
from funcAnnot.kegg.koenrich.keggenrich import enrich as koenrich
from progsuit import Configuration, Group_data, getAbsPath, matchpath
from get_gene_length import len_for_Rsp
from get_geneids import get_geneids
	
FILEPATH=os.path.realpath(__file__)
FILEDIR=os.path.dirname(FILEPATH)
CPTDE="cptDE.py"
SEQ2EXP="seq2exp.py"
FUNC="func.py"

def run_subp(argv,program):
	#program shoule be in {CPTDE,SEQ2EXP,FUNC}
	del argv[0]
	argv[0] = FILEDIR+"/"+program
	subprocess.call(argv)

if __name__ == '__main__':

	'''
	import argparse, sys
	
	parser = argparse.ArgumentParser(description='RNA-seq analyse pip')
	parser.add_argument('program',help='all for whole pip/ali for alignment/ass for assembly/\
	cptDE for compute different expression',choices=['all','cptDE','seq2exp','func'])
	parser.add_argument('-c','--conf',help='configuration file',nargs='?',\
	type=argparse.FileType('r'),default='configuration.txt')
	parser.add_argument('-1','--fq1',help='fq_1',nargs='*')
	parser.add_argument('-2','--fq2',help='fq_2',nargs='*')
	parser.add_argument('-g','--group_data',help='group_data file. conflict with -1 -2'\
	,nargs='?',type=argparse.FileType('r'))
	parser.add_argument('-s','--samples_table',help='cuffnorm samples.table file'\
	,nargs='?')
	parser.add_argument('-e','--expressionf',help='when chosing cptDE, this option is required.'\
	,nargs='?')
	parser.add_argument('-d','--DEfile',help='DEfile(NOISeq out). when chosing func,\
	this option is required.',nargs='?')
	parser.add_argument('-t','--template_gff',help='when chosing cptDE or func,\
	this option is used.',nargs='?')
	parser.add_argument('-o','--outpath',help='outpath',nargs='?')
	args = parser.parse_args(sys.argv[1:])

#	samples_table = None
#	if args.samples_table:
#		samples_table = args.samples_table
#	myconf = Configuration(args.conf)
#	ali_path = getAbsPath(args.outpath)#home path for alignment results
#	ali_name = 'mapped.sam'#alignment result name
	
#	if args.group_data and args.fq1:
#		sys.stderr.write("You can't specify both -g and -1 the same time")
#		raise StandardError
		
#	if args.group_data:
#		mygroup_data = Group_data(args.group_data)
#		fqnames, myfq1, myfq2 = mygroup_data.get2pip()
#	else:
#		myfq1 = args.fq1
#		myfq2 = args.fq2
#		fqnames = ["" for i in range(len(myfq1))]
		
#	if args.fq2 and len(args.fq1) != len(args.fq2):
#		sys.stderr.write("-2 should be as long as -1\n")
#		sys.exit(2)
	
	if args.program == 'all':
		
		gene_fpkm, htcount, quants, merged = \
		seq2exp(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data)
		nullexp,upexp,downexp = \
		cptDE(myconf,ali_path,mygroup_data,htcount,None,merged)
		func_annot(myconf,merged)
	'''	
	import argparse, sys
	
	parser = argparse.ArgumentParser(description='RNA-seq analyse pip')
	parser.add_argument('program',help='all for whole pip/ali for alignment/ass for assembly/\
	cptDE for compute different expression',choices=['all','cptDE','seq2exp','func'])
	args=parser.parse_args(sys.argv[1:2])

	if args.program == 'cptDE':
		run_subp(sys.argv,CPTDE)		

	if args.program == 'seq2exp':
		run_subp(sys.argv,SEQ2EXP)		
		
	if args.program == 'func':
		run_subp(sys.argv,FUNC)		
		
