#!/usr/bin/env python

'''
[Usage]./annot_from_db.py idfile dbfile outfile
idfile:
XLOC_00001
XLOC_00012
XLOC_00344	99	0.005
......
#only the first column be delivered to this script

dbfile:
XLOC_00001	ko00030	ko00035
XLOC_00003	ko00055
XLOC_00003	ko01004
XLOC_00003	ko00345	ko00445	ko00556
...

outfile:
XLOC_00001	ko00030
XLOC_00001	ko00030
XLOC_00334	ko00479
...
'''

def read_db_f(db_f,has_header=False,sep="\t"):

	def parse_line(line,sep,db):
		tmp = line.rstrip("\n").split(sep)
		id = tmp[0]
		tmp = tmp[1:]#del id
		db.setdefault(id,[]).extend(tmp)
	try:
		db_f = open(db_f)
	except TypeError:
		pass
	
	db = {}
	#db = {gene_id:[annot1,annot2,annot3,...],...}
	
	#parse db file

	with db_f:
		if has_header:
			hd_line = db_f.readline()
		for eachline in db_f:
			parse_line(eachline,sep,db)
	
	return db
	
def output(idfile,db,outfile,sep="\t",idindx=0,has_header=True):
	#db = {gene_id:[annot1,annot2,annot3,...],...}
	
	try:
		idfile = open(idfile)
	except TypeError:
		pass
	try:
		outfile = open(outfile,'w')
	except TypeError:
		pass
	
	with idfile:
		if has_header:
			hd_line = idfile.readline()
		for each in idfile:
			tmp = each.rstrip("\n").split(sep)
			id = tmp[idindx]
			annots = db.get(id)
			if not annots:
				continue
			for annot in annots:
				outfile.write(id+sep+annot+"\n")

def annot(idfile,dbfile,outfile,sep="\t",db_header=False,idf_header=True):

	annot_db = read_db_f(dbfile,db_header,sep)#annot db file hasn't header line
	output(idfile,annot_db,outfile,sep,0,idf_header)
	
def main(argv):

	import argparse, sys, os
	
	parser = argparse.ArgumentParser(description='annot_from_db')
	parser.add_argument('idfile',help='differential expression annalys out file,\
	typically: NOISeq out file(the first column is gene_id and has header line)',nargs='?')
	parser.add_argument('dbfile',help='go or kegg annoted db',nargs='?')
	parser.add_argument('outfile',help='annoted file',nargs='?')
	parser.add_argument('-n','--has_header',help='if idfile hasn\'t header, set this option',\
	default=True,action='store_false')
	parser.add_argument('--sep',help='sep for id file and output',default="\t")
	args = parser.parse_args(argv)
	
	annot_db = read_db_f(args.dbfile,False,args.sep)#annot db file hasn't header line
	output(args.idfile,annot_db,args.outfile,args.sep,0,args.has_header)
	
if __name__ == '__main__':
	
	import sys
	
	main(sys.argv[1:])