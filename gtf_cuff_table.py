#!/usr/bin/env python

'''
parser of gtf and cuffnorm out attribute table
_version = 0.0.0
'''

import re, sys
from collections import OrderedDict

attr_pattern = re.compile(r'\s*([^\s=;]+)[\s=]+([^"][^;]*|"[^"]*")')

def Parse_attr(pattern,attr):

    '''
    Parses Gff attribution string and results it as a dictionary
    '''

    attr_dic = {}
    tmplst = pattern.findall(attr)#there is a ";" at last each attr contain 3 element
    for key, value in tmplst:
        attr_dic[key] = value.strip("\"")
        
    return attr_dic

class Rec_Gff(object):

    '''
    one line of a standard gff file.
    '''
    
    def __init__(self,lst,chr=0,sour=1,type=2,start=3,end=4):
    
        self.chr = lst[chr]
        self.source = lst[sour]    
        self.type = lst[type]
        self.start = lst[start]
        self.end = lst[end]
        self.score = lst[5]
        self.strand = lst[6]
        self.phase = lst[7]
        self.attr = Parse_attr(pattern,lst[8])
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

class Records(object):

    _attr_pattern = re.compile(r'\s*([^\s=;]+)[\s=]+([^"][^;]*|"[^"]*")')
    
    def __init__(self,lst,source='gff'):
        #source : gff/cuff
        
        self.ID = ""
        self.parent = ""
        self.child = []
        self.annot = ""
        self.KO = ""
        self.GO = ""
        self.locus = ""
        if source == 'cuff':
            self.cuff_style(lst)
        if source == 'gff':
            self.gff_style(lst)
        
    def cuff_style(self,lst):
    
        self.ID = lst[0]
        self.class_code = lst[1]
        self.nearest_ref_id = lst[2] 
        self.gene_id = lst[3]
        self.gene_short_name = lst[4] 
        self.tss_id = lst[5]
        self.locus = lst[6]  
        self.length = lst[7]
        
    def gff_style(self,lst):
    
        self.chr = lst[0]
        self.source = lst[1]    
        self.type = lst[2]
        self.start = int(lst[3])
        self.end = int(lst[4])
        self.score = lst[5]
        self.strand = lst[6]
        self.phase = lst[7]
        self.attr = Parse_attr(Records._attr_pattern,lst[8])
        
    def get_attr(self,key):
    
        try:
            return self.attr.get(key)
        except:#need fix
            return None
            
    def write_out(self,outfile,sep="\t"):
    
        outfile.write(sep.join([self.chr,self.start,self.end,self.strand,self.attr.get("Name")]))
        
class Gene(Records):

    def __init__(self,lst,source='cuff',l_id='gene_id',l_name='Name',l_tid='transcript_id'):
        #source : gff/cuff
        #l_id/l_name/l_tid = lable of id/name/transcript id
        
        super(Gene,self).__init__(lst,source)
        if source == 'cuff':
            self.old_names = self.get_short_name(self.gene_short_name)
            self.old_name = self.get_uniq_name(self.old_names)#old gene name
            self.isf_child = []#isoform children
            self.tss_child = []#tss children
        if source == 'gff':
            self.ID = self.attr.get(l_id)
            self.name = self.attr.get(l_name)
            
    def get_short_name(self,short_name,sep=',',null='-'):
        #convert gene_short_name to list
        
        if short_name == null:
            return None
        return short_name.split(sep)
        
    def get_uniq_name(self,names):
    
        if not names:
            return None
        assert isinstance(names,list)
        if len(names) == 1:
            return ''.join(names)
        else:
            return None
            
class Isoform(Records):

    def __init__(self,lst,source='cuff'):
    
        super(Isoform,self).__init__(lst,source)
        self.parent = self.gene_id

class Tss(Records):

    def __init__(self,lst,source='cuff'):
    
        super(Tss,self).__init__(lst,source)
        self.parent = self.gene_id
        self.isf_child = []#isoform children

class Cuff_table(object):

    '''
    tracking_id     class_code      nearest_ref_id  gene_id gene_short_name tss_id  locus   length
    XLOC_000001     -       -       XLOC_000001     -       TSS1    CI01000000:32133-102423 -
    XLOC_000002     -       -       XLOC_000002     -       TSS2    CI01000000:32133-102423 -
    XLOC_000003     -       -       XLOC_000003     everse_transcriptase_[Danio_rerio]      TSS3    CI01000000:107961-111210        -
    XLOC_000004     -       -       XLOC_000004     _TPAxp:_polyprotein_[Danio_rerio]       TSS4    CI01000000:111528-115356        -
    '''
    def __init__(self,mytable,type='gene',sep="\t",has_header=True):
        #type = 'gene'/ 'isoform'/ 'tss'
        #mytable = cuff_table file object
        
        self._data = OrderedDict()
        self.header = ''
        self.type = self.loadtype(type)
        self.loadfile(mytable,sep,has_header)
        
    def loadfile(self,infile,sep="\t",has_header=True):
    
        if has_header:
            self.header = infile.readline()
        for eachline in infile:
            tmp = eachline.strip().split(sep)
            self._data[tmp[0]] = self.load_rec(tmp)
            
    def loadtype(self,type):
    
        if type == 'gene' or type == 'isoform' or type == 'tss':
            return type
        return None
    
    def load_rec(self,lst):
    
        if not self.type:
            return Records(lst,'cuff')
        if self.type == 'gene':
            return Gene(lst)
        if self.type == 'isoform':
            return Isoform(lst)
        if self.type == 'tss':
            return Tss(lst)
    
    def get(self,ID):

        return self._data.get(ID)
        
class Gff(object):

    '''
    convert a gff file to an iterator.
    and provide gene id trace method (get_gene()) when make_store set True
    '''

    def __init__(self,file,sep="\t",has_header=False,flag='cuff',make_store=False,\
    selection=None,selections=[]):
        #flag: change the lable of gene_id attribution
        #selection(s): chose one type(s) to be stored or iteration.
        
        if selection:
            self.selections = set([selection])
        elif selections:
            self.selections = set(selections)
        else:
            self.selections = None
        
        self.has_header = has_header
        self.sep = sep
        self.l_id = 'ID'
        #l_id lable of id
        if flag == 'cuff':
            self.l_id = 'gene_id'
        
        self.make_store = make_store    
        if make_store:
            self._data = []
            self._genes = OrderedDict()
            self._names = OrderedDict()
            self._loadfile(file)
        else:
            self.infile = file
        
    def __iter__(self):
    
        if self.make_store:
            return self._next()
        else:
            return self._next_unstored()
    
    def __contains__(self,key,flag='gene_id'):
    
        if flag == 'gene_id':
            return key in self._genes
    
    def _next(self):
    
        for x in self._data:
            yield x
            
    def _next_unstored(self):
    
        if self.has_header:
            self.header = self.infile.readline()
        for eachline in self.infile:
            tmp_rec = self._loadline(eachline)
            if tmp_rec:
                yield tmp_rec
            
    def _loadfile(self,file):
    
        if self.has_header:
            self.header = infile.readline().strip()
        indx = 0
        for eachline in file:
            tmp = self._loadline(eachline)
            if tmp:
                self._data.append(tmp)
                self._append_gene(tmp,indx)
            #self.append_name(tmp,indx)
                indx += 1
            
    def _loadline(self,line):
        
        rec = Records(line.strip().split(self.sep),'gff')
        if not self.selections:
            return rec
        elif rec.type in self.selections:
            return rec
        else:
            return None
        
    def _append_gene(self,myrec,indx):
        
        if myrec.type != "gene":
            return None
        self._append_attr(myrec,indx,self.l_id,self._genes)
            
    def append_name(self,myrec,indx):

        if myrec.type != "gene":
            return None
        self._append_attr(myrec,indx,"Name",self.names)
        
    def _append_attr(self,myrec,indx,attr_name,mydic):
        #myrec is a Records
        #indx : index in self._data
        #attr : which attribute to be the key
        #mydic : an index dict which map attr to index in _data
        
        myattr = myrec.attr.get(attr_name)
        if not myattr:
            return None
        if myattr in mydic:
            sys.stderr.write("warning: attribution %s not unique\n"%myattr)
            self.get_gene(myattr,'Name').write_out(sys.stderr)
            myrec.write_out(sys.stderr)
        mydic[myattr] = indx
    
    def _get_attr(self,key,mydic):
        #should be renamed
    
        indx = mydic.get(key)
        if indx == 0:
            return self._data[indx]
        if not indx:
            return None
        return self._data[indx]
    
    def get_gene(self,key,flag='ID'):    
        #flag: ID/Name
        
        if not self.make_store:
            sys.stderr.write("get_gene can only be used when 'make_store' set as True")
            raise EOFError
            
        if flag == 'ID':
            return self._get_attr(key,self._genes)
        if flag == 'Name':#can not be used    
            return self._get_attr(key,self._names)
            
    def get_genes(self):
    
        if not self.make_store:
            sys.stderr.write("get_genes can only be used when 'make_store' set as True")
            raise EOFError
            
        return self._genes.values()
        
def main():
    
    import sys

    gc_final = open(sys.argv[1])
    merged = open(sys.argv[2])
    guidef = open(sys.argv[3])
    #mygene = Cuff_table(fp)
    #print(mygene.get("XLOC_000005"))
    #print(mygene.type)
    myguide = Gff(guidef,flag='cuff')
    mymerge = Gff(merged,flag='cuff')
    mygc = Gff(gc_final,selection='gene',make_store=True,flag='gff')
    
    curr_id = ''
    coord = [0,0]
    mymerge.gene_edge = {}

    def update_edge(crd1,crd2):
        if not crd2:
            newcrd = crd1
        else:
            newcrd = [min(crd1[0],crd2[0]),max(crd1[1],crd2[1])]
        return newcrd
        
    for each in mymerge:#make mymerge.gene_edge dictionary
        gene_id = each.get_attr("gene_id")
        if gene_id != curr_id:
            if not curr_id:
                curr_id = gene_id
                coord = [each.start,each.end]
                continue
            old_edge = mymerge.gene_edge.get(curr_id)
            mymerge.gene_edge[curr_id] = update_edge(coord,old_edge)
            curr_id = gene_id
            coord = [each.start,each.end]
            continue
        #gene_id == curr_id:
        coord = update_edge(coord,[each.start,each.end])
    old_edge = mymerge.gene_edge.get(curr_id)
    mymerge.gene_edge[curr_id] = update_edge(coord,old_edge)
            
    for each in myguide:
        #print(each.get_attr("gene_id"))
        gene_id = each.get_attr("gene_id")
        old_id = each.get_attr("oId")
        oldgene = mygc.get_gene(old_id)
        estart,eend = mymerge.gene_edge.get(gene_id)#edge of start/end
        d = float(min(eend,each.end)-max(estart,each.start))
        p1 = d/(each.end-each.start)
        p2 = d/(eend-estart)
        sys.stdout.write("\t".join(map(str,[gene_id,old_id,each.start,each.end,\
        oldgene.start,oldgene.end,estart,eend,p1,p2]))+"\n")
        #len_diff = int(oldgene.end) - int(oldgene.start) - int(each.end) + int(each.start)
        #sys.stdout.write(str(len_diff))
        
if __name__ == '__main__':

    import sys
    sys.exit(main())
