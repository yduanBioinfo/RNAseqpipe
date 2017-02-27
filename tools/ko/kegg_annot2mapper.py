#!/usr/bin/env python

import sys
'''
this script is used to generate a file which format is fit for keggmapper.
source file will be a kegg annotation file and a file including geneid and logFC.

example: ./kegg_annot2mapper.py annotfile expfile outfile
'''

annotfile = open(sys.argv[1])
expfile = open(sys.argv[2])
outfile = open(sys.argv[3],"w")

expf_ID = 0; expf_logFC = 5#the position of geneid and logFC
expf_has_header = True; annotf_has_header = False
if expf_has_header:
	expfile.readline()
if annotf_has_header:
	annotfile.readline()

ko_expdata = {}; expdata = {}
#ko_expdata = {ko1:value,ko2:value,ko3:value,...}
#the value can be sum of logFC or sum of sign
#expdata = {ID1:logFC,ID2:logFC,ID3:logFC,...}

def get_expdata(file,ID,FC,sep="\t"):
	
	'''ID or FC represent which column is the ID or logFC to be
	'''
	data = {}#expdata={ID1:logFC,ID2:logFC,ID3:logFC,...}
	for eachline in file:
		tmp = eachline.rstrip().split(sep)
		if tmp[ID] in data:
			print("Warning:%s is not unique,only the first is counted"% tmp[ID])
			continue
		data[tmp[ID]] = float(tmp[FC])
	return data


def get_ko_exp(annotfile,expdata,method="logFC",sep="\t"):

	'''method can be "logFC" or "sign"
	annotfile:
	Gene1	KO1
	Gene2	KO2
	Gene3	KO3
	'''
	data = {}#ko_expdata={ko1:value,ko2:value,ko3:value,...}
	for eachline in annotfile:
		tmp = eachline.rstrip().split(sep)
		if len(tmp) != 2:
			continue
		myexp = expdata[tmp[0]]
		if method == "sign":
			if myexp < 0:
				myexp = -1
			elif myexp > 0:
				myexp = 1
		data.setdefault(tmp[1],0)
		data[tmp[1]] += myexp
	return data

expdata=get_expdata(expfile,expf_ID,expf_logFC)
ko_expdata=get_ko_exp(annotfile,expdata,'sign')

def write_mapperfile(data,outfile,upbColor='red',upfColor='black',\
downbColor='green',downfColor='black',unbColor='blue',unfColor='black'):

	'''
	data = {KO1:exp,KO2:exp,KO3:exp,...}
	'''
	for key, value in data.items():
		outfile.write(key+"\t")
		if value > 0:
			outfile.write("%s,%s\n"%(upbColor,upfColor))
		elif value < 0:
			outfile.write("%s,%s\n"%(downbColor,downfColor))
		else:
			outfile.write("%s,%s\n"%(unbColor,unfColor))
		
write_mapperfile(ko_expdata,outfile)