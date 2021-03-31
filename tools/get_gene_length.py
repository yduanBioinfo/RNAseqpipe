#!/usr/bin/env python

#from RNAseqpipe.gtf_cuff_table import Gff
from collections import OrderedDict as Ordic
from dypylib.bio.seq.base import Gff

'''
input gtf should only contain exon feature.
'''

def make_feature_dict(mygff,key="transcript_id"):

	data={}
	for rec in mygff:
		data.setdefault(rec.attr.get(key),[]).append(rec)
	return data

def stat_db(mydb):
	#be fit for transcripts survey
	#data:(transcript name, exon counts, transcript length)
	#exon_length:length of each exon

	data=[]
	exon_length = []
	for key,value in mydb.items():
		count=len(value)
		exons_l=map(len,value)#length of each exon
		exon_length.extend(exons_l)
		length=sum(exons_l)
		data.append((key,count,length))
	return data,exon_length

class Gff_l(Gff):

    '''
    add get length function to Gff class 
    
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
    '''
#    def attr_pos(self,t_attr='gene_id'):

 #       pass

    def make_feature_dict(self,key="transcript_id",selected=["exon"]):

        data=Ordic()
        for rec in self:
            if rec.type not in selected:
                continue
            data.setdefault(rec.attr.get(key),[]).append(rec)
        return data

    def attr_lener(self,t_attr='gene_id'):
    #return iterator of (id,length)

        tdata = self.make_feature_dict(t_attr)
        
        for key,value in tdata.items():
            exons_l=map(len,value)#length of each exon
            length=sum(exons_l)
            yield key,length

#    def attr_len(self,t_attr='gene_id'):
        #return iterator of (id,length)
        
#        outdic = Ordic()
#        for key, value in self.attr_pos(t_attr).items():
#            outdic[key] = value[1] - value[0]
#        return outdic
    
    def attr_len(self,t_attr='gene_id'):
    #return data of {id,length}
        
        data = Ordic()
        for key, value in attr_lener(t_attr):
            data[key]=value
        return data

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
    parser.add_argument('-t', '--t_attr', help='id attribute name',nargs='?',default='gene_id')
    parser.add_argument('-o','--outfile',nargs='?',help='outfile default: stdout',\
    default=sys.stdout,type=argparse.FileType('w'))
    args = parser.parse_args(argv)

    #mygff = Gff_l(args.gff)
    mygff = Gff_l(args.gff)
    write_length_f(mygff.attr_lener(args.t_attr),args.outfile)
        
if __name__ == '__main__':

    import sys
    
    main(sys.argv[1:])
