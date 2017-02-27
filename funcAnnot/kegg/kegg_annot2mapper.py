#!/usr/bin/env python

import sys, math
'''
this script is used to generate a file which format is fit for keggmapper.
source file will be a kegg annotation file and a file including geneid and logFC.
_version = 0.1.0
example: ./kegg_annot2mapper.py annotfile expfile outfile
'''

#annotfile = open(sys.argv[1])
#expfile = open(sys.argv[2])
#outfile = open(sys.argv[3],"w")

#expf_ID = 0; expf_logFC = 5#the position of geneid and logFC
#expf_has_header = True; annotf_has_header = False
#if expf_has_header:
#	expfile.readline()
#if annotf_has_header:
#	annotfile.readline()

#ko_expdata = {}; expdata = {}
#ko_expdata = {ko1:value,ko2:value,ko3:value,...}
#the value can be sum of logFC or sum of sign
#expdata = {ID1:logFC,ID2:logFC,ID3:logFC,...}

def get_expdata(file,ID,FC,sep="\t"):
	
	'''ID or FC represent which column is the ID or logFC to be
	'''
	data = {}#expdata={ID1:logFC,ID2:logFC,ID3:logFC,...}
	for eachline in file:
		tmp = eachline.rstrip().split(sep)
		myid = tmp[ID].strip()
		if myid in data:
			print("Warning:%s is not unique,only the first is counted"% tmp[ID])
			continue
		data[myid] = float(tmp[FC])
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

#expdata=get_expdata(expfile,expf_ID,expf_logFC)
#ko_expdata=get_ko_exp(annotfile,expdata,'sign')

def cmpt_zsc(annotfile,expdata,sep="\t"):#compute_zscores

	'''
	zscores = (up-down)/squr(n)
	up/down represent the count of up/down regulation genes
	n = up + down

	annotfile:
	Gene1	KO1
	Gene2	KO2
	Gene3	KO3
	'''
	data = {}#ko_expdata={ko1:[up_count,down_count],ko2:[up,down],ko3:[up,down],...}
	for eachline in annotfile:
		tmp = eachline.rstrip().split(sep)
		if len(tmp) != 2:
			continue
		myexp = expdata.get(tmp[0])
		if not myexp:
			continue
		data.setdefault(tmp[1],[0,0])
		if myexp > 0:
			data[tmp[1]][0] += 1
		if myexp < 0:
			data[tmp[1]][1] += 1
	outdata = {}
	for key, value in data.items():
		if value[0] + value[1] == 0:
			continue
		outdata[key] = (value[0] - value[1])/math.sqrt(value[0]+value[1])
	return outdata

def write_mapperfile(data,outfile,upbColor='red',upfColor='black',\
downbColor='green',downfColor='black',unbColor='blue',unfColor='black'):

	'''
	data = {KO1:exp,KO2:exp,KO3:exp,...}
	'''
	for key, value in data.items():
		outfile.write(key+"\t")
		if value > 0.5:
			print(value)
			outfile.write("%s,%s\n"%(upbColor,upfColor))
		elif value < -0.5:
			print(value)
			outfile.write("%s,%s\n"%(downbColor,downfColor))
		else:
			outfile.write("%s,%s\n"%(unbColor,unfColor))

#ko_expdata=cmpt_zsc(annotfile,expdata)
#write_mapperfile(ko_expdata,outfile)

def main(argv):

	import argparse
	parser = argparse.ArgumentParser(description='keggmapper suit.\
	this script is used to generate a file which format is fit for keggmapper.\
	source file will be a kegg annotation file and a file including geneid and logFC.')
	parser.add_argument('annotfile',help='',nargs='?',type=argparse.FileType('r'))
	parser.add_argument('expfile',help='col0:ID;col5:logFC',nargs='?',type=argparse.FileType('r'))
	parser.add_argument('outfile',help='',nargs='?',type=argparse.FileType('w'))
	args = parser.parse_args(argv)

	expf_ID = 0; expf_logFC = 5#the position of geneid and logFC
	expf_has_header = True; annotf_has_header = False
	if expf_has_header:
		args.expfile.readline()
	if annotf_has_header:
		args.annotfile.readline()

	ko_expdata = {}; expdata = {}
	#ko_expdata = {ko1:value,ko2:value,ko3:value,...}
	#the value can be sum of logFC or sum of sign
	#expdata = {ID1:logFC,ID2:logFC,ID3:logFC,...}

	expdata=get_expdata(args.expfile,expf_ID,expf_logFC)
	#ko_expdata=get_ko_exp(annotfile,expdata,'sign')
	
	ko_expdata=cmpt_zsc(args.annotfile,expdata)
	write_mapperfile(ko_expdata,args.outfile)

if __name__ == '__main__':

	main(sys.argv[1:])

