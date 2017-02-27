#!/usr/bin/env python

from gtf_cuff_table import Gff
from collections import OrderedDict as Ordic

class Gff_l(Gff):

	'''
	add get length function to Gff class 
	'''
	
	def attr_pos(self,t_attr='gene_id'):
		#t_attr: type of attr.
		#return a dict of attr_id:[min_start,max_end]
		
		if self.make_store == False:
		
			outdic = Ordic()
			for each_rec in self:
				id = each_rec.attr[t_attr]
				st = int(each_rec.start)
				ed = int(each_rec.end)
				oldpos = outdic.get(id)
				if oldpos:#update position list
					st = min(oldpos[0],st)
					ed = max(oldpos[1],ed)
				outdic[id] = [st,ed]
			return outdic
			
		else:
			sys.stderr.write("stored Gff don\'t support this function till now.")
			raise TypeError

	def attr_len(self,t_attr='gene_id'):
		#return iterator of (id,length)
		
		outdic = Ordic()
		for key, value in self.attr_pos(t_attr).items():
			outdic[key] = value[1] - value[0]
		return outdic
	
	def attr_lener(self,t_attr='gene_id'):
		#return iterator of (id,length)
		
		for key, value in self.attr_pos(t_attr).items():
			yield key, value[1]-value[0]

	def get_length_array(self,t_attr='gene_id'):
		
		outlst = []
		for key, value in self.attr_lener(t_attr):
			outlst.append((key,value))
		return outlst
		
def write_length_f(length_dic,outfile,sep='\t'):
	#length_dic: is actully a iterator of (id, length) pair.

	for id, length in length_dic:
		outfile.write(id+sep+str(length)+"\n")

def len_for_Rsp(gff,t_attr='gene_id'):
	#get length for RNAseqpip
	mygff = Gff_l(gff)
	return mygff.get_length_array(t_attr)
	
def main(argv):

	import argparse, sys
	
	parser = argparse.ArgumentParser(description='merge htseq-count out file')
	parser.add_argument('gff',help='gff files',nargs='?',type=argparse.FileType('r'))
	parser.add_argument('-o','--outfile',nargs='?',help='outfile default: stdout',\
	default=sys.stdout,type=argparse.FileType('w'))
	args = parser.parse_args(argv)

	mygff = Gff_l(args.gff)
	write_length_f(mygff.attr_lener(),args.outfile)
		
if __name__ == '__main__':

	import sys
	
	main(sys.argv[1:])