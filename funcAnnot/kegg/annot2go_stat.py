#!/usr/bin/env python

import sys

'''This script convert go annot file to the file format go_stats required 
'''

def annot2go_stat(annotfile,stats_file,symble="ISA"):
	#symble = "ISA"#accoding to Gene ontology
	
	annotfile = open(annotfile)
	stats_file = open(stats_file,'w')
	
	for each in annotfile:
		tmp = each.strip().split("\t")
		stats_file.write("%s\t%s\t%s\n"%(tmp[1],symble,tmp[0]))
	
	annotfile.close()
	stats_file.close()
	
if __name__ == '__main__':

	import sys
	
	annotfile = sys.argv[1]
	stats_file = sys.argv[2]
	annot2go_stat(annotfile,stats_file)