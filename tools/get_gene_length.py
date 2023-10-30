#!/usr/bin/env python

#from RNAseqpipe.gtf_cuff_table import Gff
from collections import OrderedDict as Ordic
from dypylib.bio.seq.base import Gff
from dypylib.bio.seq.Annotation import Genome

''' Input should be GTF.
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
    ''' Class of collections of methods for getting length genome features in GFF/GTF file.
    Input should be a file.
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

def get_length_array(myGTF, transcripts=False, feature_id='gene_id'):
    """ Get list of [Gene/Tx ID, length]
    myGTF is GtfDict 
    transcripts: True, get Transcripts information, otherwise Gene.
    feature_id: Which feature_id is used as key. Important: it can output 
        gene_names/gene_id/gene_symbol with gene level length.
    """
    outdata = []
    if transcripts:
        dict_txs = myGTF.get_TxDict()
        for tx in dict_txs.values():
            outdata.append([tx.ID,len(tx)])
    else:
        dict_genes = myGTF.get_GeneDict()
        for gene in dict_genes.values():
            outdata.append([gene.get_attribute(feature_id),len(get_represent_tx(gene))])
    return outdata

def _get_length_array(gff, t_attr):
    """ Enable len_for_Rsp function run with new engine """
    
    transcripts = False
    if t_attr in ("tx_id", "transcript_id"):
        transcripts = True
    return get_length_array(gff, transcripts)

def len_for_Rsp(gff,t_attr='gene_id'):
    """ get length for RNAseqpip"""
    myGTF = Genome(args.gff)
    return _get_length_array(myGTF, t_attr)

def len_for_Rsp_old(gff,t_attr='gene_id'):
    #get length for RNAseqpip abord
    mygff = Gff_l(gff)
    return mygff.get_length_array(t_attr)

def get_represent_tx(gene):
    """ Get representative transcript of one gene, the longest one.
        This function should be incorporated into dypylib.
    """
    max_len = 0
    for tx in gene.values():
        if len(tx) > max_len:
            rep_tx = tx
            max_len = len(tx)
    return rep_tx
    
def main(argv):

    import argparse, sys
    
    parser = argparse.ArgumentParser(description='Get sequence length of each gene/transcript')
    parser.add_argument('gff',help='GTF file (GFF file is not supported)',nargs='?',type=argparse.FileType('r'))
    parser.add_argument('-t', '--t_attr', help='Attribute name for output',nargs='?',default='gene_id')
    parser.add_argument('--transcripts', help='Get length of transcript (default: gene)', action='store_true')
    parser.add_argument('-o','--outfile',nargs='?',help='outfile default: stdout',\
    default=sys.stdout,type=argparse.FileType('w'))
    args = parser.parse_args(argv)

    myGTF = Genome(args.gff)
    for ID, length in get_length_array(myGTF, args.transcripts,args.t_attr):
        args.outfile.write("{}\t{}\n".format(ID, length))
        
if __name__ == '__main__':

    import sys
    
    main(sys.argv[1:])
