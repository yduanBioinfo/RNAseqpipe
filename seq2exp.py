#!/usr/bin/env python

from __future__ import print_function
import sys

from RNAseqpipe import hisat, cufflinks, htseq, verse, stringtie
from RNAseqpipe.expression.salmon_quantify import salmonpip
from RNAseqpipe.progsuit import Configuration, Group_data, getAbsPath, matchpath
from RNAseqpipe.run_RNAseqpipe import log, add_arguments, BASE_CONF

def seq2exp(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data,run_cdiff=True):
    pipe = myconf["all"].get("pipe")
    if pipe == "hch":
        # Hisat+ Cufflinks+ HTseq
        return hisat_cufflinks_htseq(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data,run_cdiff)

    if pipe == "hcv":
        return hisat_cufflinks_verse(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data,run_cdiff)

    if pipe == "hsh":
        # Hisat + Stringtie + HTseq
        return hisat_stringtie_htseq(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data)

    if pipe == "hsv":
        # Hisat + Stringtie + Verse
        return hisat_stringtie_verse(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data)

    if pipe == "hh":
        # Hisat + Htseq
        return hisat_htseq(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data)
    if pipe == "hv":
        # Hisat + Verse
        return hisat_verse(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data)

def sub_hisat_cufflinks(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data,run_cdiff=True):
    # sub_* sub processing should be called repeatly.
    # from fastq to assembly gff
    # run_cdiff: if cuffdiff step should be performed

    ali_ress,sort_ress = hisat.pip_hisats(myconf,myfq1,myfq2,fqnames,ali_path,ali_name)
    #alignment results,samtools sorted results

    cufflinks_ress = cufflinks.cufflinks(myconf,sort_ress)#cufflinks results

    #write path/assembly.txt
        
    assembly = ali_path+"/assembly.txt"
    cu_ress_f = open(assembly,'w')
    cu_ress_f.write("\n".join(cufflinks_ress))
    cu_ress_f.close()
        
    merged = cufflinks.cuffmerge(myconf,assembly,ali_path+"/merged_asm")
    quants = cufflinks.cuffquants(myconf,sort_ress,merged)
    norm = cufflinks.cuffnorm(myconf,quants,merged,ali_path)
    gene_fpkm = norm+"/genes.fpkm_table"#cufflinks fpkm
    samples_table = norm+"/samples.table"
    #run cuffdiff
    treat_group = mygroup_data.get_g_cdgp()
    #{group1:[dataname1,dataname2,...],group2:[],...}
    try:
        gpnames1,gpnames2=treat_group.values()
    except:
        print(treat_group.values(),file=sys.stderr)
        run_cdiff = False
        log.error("Warning:only two group data can run cuffdiff!")
        log.error("cuffdiff step skipped")
    
    if run_cdiff:
        quants1,quants2 = matchpath(gpnames1,gpnames2,quants)
        cu_diff = ali_path+"/diffout"
        cufflinks.cuffdiff(myconf,merged,quants1,quants2,cu_diff)
    
    #cuffdiff end
    #return expressionf    
    return gene_fpkm, sort_ress, quants, merged

def hisat_cufflinks_htseq(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data,run_cdiff=True):
    #from fastq to assembly gff
    #run_cdiff: if cuffdiff step should be performed

    gene_fpkm, sort_ress, quants, merged = sub_hisat_cufflinks(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data,run_cdiff)
    #cuffdiff end
    htcount = htseq.htseqpip(myconf,sort_ress,ali_path,merged,plotqa=False)
    #return expressionf    
    return gene_fpkm, htcount, quants, merged

def hisat_stringtie_htseq(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data):
    ali_ress,sort_ress = hisat.pip_hisats(myconf,myfq1,myfq2,fqnames,ali_path,ali_name)
    stringtie_ress = stringtie.stringties(myconf,sort_ress)
    # stringtie results, gtf files
    merged = stringtie.merge(myconf,stringtie_ress,ali_path+"/merged_asm/merged.gtf")
    htcount = htseq.htseqpip(myconf,sort_ress,ali_path,merged,plotqa=False)
    return htcount, merged

def hisat_htseq(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data):
    #from fastq to assembly gff
    #Don't finding novel gene.
    #Only get expression.

    ali_ress,sort_ress = hisat.pip_hisats(myconf,myfq1,myfq2,fqnames,ali_path,ali_name)
    #Target gff
    merged = myconf["all"].get("gff")

    htcount = htseq.htseqpip(myconf,sort_ress,ali_path,merged,plotqa=False)
    return htcount, merged

def hisat_cufflinks_verse(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data,run_cdiff):
    ali_ress,sort_ress = hisat.pip_hisats(myconf,myfq1,myfq2,fqnames,ali_path,ali_name)
    stringtie_ress = stringtie.stringties(myconf,sort_ress)
    # stringtie results, gtf files
    merged = stringtie.merge(myconf,stringtie_ress,ali_path+"/merged_asm/merged.gtf")
    vscount = verse.versepip(myconf,sort_ress,ali_path,merged)
    salmonpip(myconf,myfq1,myfq2,fqnames,ali_path+"/salmon",merged)
    return vscount, merged

def hisat_stringtie_verse(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data):
    ali_ress,sort_ress = hisat.pip_hisats(myconf,myfq1,myfq2,fqnames,ali_path,ali_name)
    stringtie_ress = stringtie.stringties(myconf,sort_ress)
    # stringtie results, gtf files
    merged = stringtie.merge(myconf,stringtie_ress,ali_path+"/merged_asm/merged.gtf")
    vscount = verse.versepip(myconf,sort_ress,ali_path,merged)
    salmonpip(myconf,myfq1,myfq2,fqnames,ali_path+"/salmon",merged)
    return vscount, merged

def hisat_verse(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data):
    #from fastq to assembly gff
    #Don't finding novel gene.
    #Only get expression.

    log.debug("Running hisat_verse")
    ali_ress,sort_ress = hisat.pip_hisats(myconf,myfq1,myfq2,fqnames,ali_path,ali_name)
    #Target gff
    vmerged = myconf["all"].get("gtf")
    if vmerged is None:
        log.warning("Verse only support gtf, which is not provied.\nTry with gff.")
        vmerged = myconf["all"].get("gff")

    smerged = myconf["all"].get("gff")
    if smerged is None:
        log.warning("Gffread are used in salmon. Gffread work better on gff,\
                which is not provied.\nTry with gtf.")
        smerged = myconf["all"].get("gtf")

    vscount = verse.versepip(myconf,sort_ress,ali_path,vmerged)
    salmonpip(myconf,myfq1,myfq2,fqnames,ali_path+"/salmon",smerged)
    return vscount, smerged

def count2exp(myconf,count,outdir=None):
    if outdir == None:
        # countdir --> outdir
        pass

def print_conf(conf):
    for k,v in conf.items():
        print("top level")
        print(k)
        for k1,v1 in v.items():
            print(v1)
    
def main(argv): 
    import argparse, sys
    
    parser = argparse.ArgumentParser(description='RNA-seq seq2exp program')
    parser.add_argument('-1','--fq1',help='fq_1',nargs='*')
    parser.add_argument('-2','--fq2',help='fq_2',nargs='*')
    add_arguments(parser)
    args = parser.parse_args(argv[1:])

    if len(argv) == 1:
        parser.print_usage()
        sys.exit(1)

    myconf = Configuration(args.conf, base_conf=BASE_CONF)
    ali_path = getAbsPath(args.outpath)#home path for alignment results
    ali_name = 'mapped.sam'#alignment result name

    #set log
    if args.verbose and args.quite:
        log.error("You can only choose one option between verbose and quite!")
        sys.exit(1)
    #waning or higher level
    log.setLevel(20)
    #output all
    if args.verbose:
        log.setLevel(1)
    #Error or higher level
    if args.quite:
        log.setLevel(40)
    
    if args.group_data and args.fq1:
        log.error("You can't specify both -g and -1 the same time")
        sys.exit(1)
       # raise StandardError("You can't specify both -g and -1 the same time")
    #input data processing
    #convert group_data 2 paired data
    #check if fq1 contains as many files as fq2 is when dealing with paired data
    if args.group_data:
        mygroup_data = Group_data(args.group_data)
        fqnames, myfq1, myfq2 = mygroup_data.get2pip()
    else:
        myfq1 = args.fq1
        myfq2 = args.fq2
        fqnames = ["" for i in range(len(myfq1))]
        
    if args.fq2 and len(args.fq1) != len(args.fq2):
        log.error("-2 should be as long as -1")
        sys.exit(1)
       # raise StandardError("-2 should be as long as -1\n")

    seq2exp(myconf,myfq1,myfq2,fqnames,ali_path,ali_name,mygroup_data,run_cdiff=False)
        
if __name__ == '__main__':

    import sys
    main(sys.argv)
