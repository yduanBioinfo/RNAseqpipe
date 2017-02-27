#!/usr/bin/env python

import re, string

'''
This script is designed for reading gene information from a standard gtf file,
and getting its sequence from a genome file.
3 kind of objects exists in this script.
__author : Duan-You
__version : 0.0.2
'''

class Fasta(object):
	
	'''
	data = {key:[[info],val]}
	key = geneID
	info = gene infomation
	val = sequence
	'''

	def __init__(self,file):

		self.__length = 0
		self.__data = {}
		self.loadfile(file)
	
	def __len__(self):

		return self.__length

	def __getitem__(self,key):

		return self.__data.get(key)

	def getSeq(self,key):

		return self.__data.get(key)[1]

	def getNames(self):

		return self.__data.keys()

	def loadfile(self,file):

		currseq = ''
		for eachline in file:
			if eachline[0] == ">":
				if currseq:
					self.__length += 1
					self.__data[currinfo[0]].append(currseq)
					currseq = ''
				currinfo = eachline.lstrip(">").split()
				self.__data.setdefault(currinfo[0],[]).append(currinfo[1:])
			else:
				currseq += eachline.strip()
		self.__data[currinfo[0]].append(currseq)
		self.__length += 1

	def getRC(self,seq):

		'''
		get reverse compliment sequence
		'''

		intab = 'AaTtGgCc'
		outtab = 'TtAaCcGg'
		transtab = string.maketrans(intab,outtab)
		return seq[::-1].translate(transtab)
	
	def getData(self):

		return self.__data

attr_partern = re.compile(r"\s*([^\s=]+)[\s=]+(.*)")

def Parse_attr(attr):

	'''Parses Gff attribution string and results it as a dictionary
	'''

	attr_dic = {}
	tmplst = map(lambda x:attr_partern.match(x),attr.strip().split(";")[:-1])#there is a ";" at last
	for each in tmplst:
		if not each:
			print("gff type error")
		val = each.group(2)
		if val.startswith('"') and val.endswith('"'):
			val = val[1:-1]
		attr_dic.setdefault(each.group(1),val)
	return attr_dic

class Gene(object):

	def __init__(self,lst,chr=0,sour=1,type=2,start=3,end=4):
	
		self.chr = lst[chr]
		self.source = lst[sour]	
		self.type = lst[type]
		self.start = lst[start]
		self.end = lst[end]
		self.score = lst[5]
		self.strand = lst[6]
		self.phase = lst[7]
		self.attr = Parse_attr(lst[8])
		self.seq = ''

	def __getattr__(self,name):

		print("%s don\'t have an attribute named %s"% (self.attr,name))

	def __len__(self):

		return int(self.end)-int(self.start)+1

	def loadseq(self,genome):

		self.seq = genome[self.chr][1][int(self.start)-1:int(self.end)]

	def getGtf(self):

		return [self.chr,self.source,self.type,self.start,self.end,\
			self.score,self.strand,self.phase,self.attr]

	def getAll(self):

		return self.getGtf.append(self.seq)

	def writeFile1(self,file):

		file.write("\t".join([">"+self.chr,self.attr,self.start,\
		self.end,self.strand])+"\n")
		file.write("\n".join(self.seq[i*70:i*70+70] for i in range(len(self.seq)/70+1))+"\n")

	def writeFile2(self,file,key):

		'''
		only write one attr
		'''
		file.write(">"+self.attr[key]+"\n")
		file.write("\n".join(self.seq[i*70:i*70+70] for i in range(len(self.seq)/70+1))+"\n")

class Gtf(object):

	'''need fix
	'''

	def __init__(self,file):

		'''
		data = [gene1,gene2,gene3,...]
		'''
		
		self.data = []
		self.loadfile(file)

	def loadfile(self,file):

		for eachline in file:
			self.data.append(Gene(eachline.strip().split("\t")))

def writeFasta(chr,seq,file):

	'''
	universe
	'''
	file.write(">"+chr+"\n")
	file.write("\n".join(seq[i*70:i*70+70] for i in range((len(seq)-1)/70+1))+"\n")

if __name__ == '__main__':

	import argparse, sys

	parser = argparse.ArgumentParser(description='find gene seq in genome')
	parser.add_argument('fastaFile',nargs='?',type=argparse.FileType('r'),\
			    default=sys.stdin,help='fastaFile [default: stdin]')
	parser.add_argument('GtfFile',nargs='?',type=argparse.FileType('r'),help='')
	parser.add_argument('IDfile',nargs='?',type=argparse.FileType('r'),help='')
	parser.add_argument('-o','--outFile',nargs='?',type=argparse.FileType('w'),\
			    default=sys.stdout,help='outfile [default: stdout]')
	args = parser.parse_args(sys.argv[1:])

	'''
	find gene seq and generate a fasta file.
	genes ID in IDfile
	get location in gtffile
	and return sequence in fastafile
	'''
	myIDs = [i.strip() for i in args.IDfile]
	mygenome = Fasta(args.fastaFile)
	curr_seq = ''
	curr_rec = ['','']
	for each in args.GtfFile:
		a = Gene(each.strip().split("\t"))
		if a.attr['gene_id'] not in myIDs:
			continue
		if a.attr['gene_id'] == curr_rec[0] and \
		a.attr['transcript_id'] != curr_rec[1]:
			continue
		if a.attr['gene_id'] != curr_rec[0] and curr_seq:
			writeFasta(curr_rec[0],curr_seq,args.outFile)
			curr_seq = ''
		curr_rec = [a.attr['gene_id'],a.attr['transcript_id']]
		a.loadseq(mygenome)
		curr_seq += a.seq
	writeFasta(curr_rec[0],curr_seq,args.outFile)
