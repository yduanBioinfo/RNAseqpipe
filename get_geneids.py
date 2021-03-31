#!/usr/bin/env python

from RNAseqpipe.gtf_cuff_table import Gff
import sys

'''usage: ./get_geneids.py gff outfile
    or ./get_geneids.py gff > outfile
    
result file :
    gene_id
XLOC_000001
XLOC_000002
XLOC_000003
XLOC_000004
XLOC_000005
XLOC_000006
XLOC_000007
XLOC_000008
XLOC_000009
...
...

'''
def get_geneids(gff,outfile):

    argvs = [gff,outfile]
    main(argvs)
    
    return outfile
    
def main(argvs):
    
    gff = argvs[0]
    try:
        gff = open(gff)
    except:
        pass
    try:
        outfile = argvs[1]
        outfilename = outfile
        outfile = open(outfile,'w')
    except:
        outfile = sys.stdout

    mygff = Gff(gff)

    outfile.write("gene_id\n")
    curr_rec = ''

    for rec in mygff:
    
        rec = rec.attr["gene_id"]
        if curr_rec == rec:
            continue
        
        outfile.write(rec+"\n")
        curr_rec = rec
    
if __name__ == '__main__':

    import sys
    
    get_geneids(sys.argv[1:])
