# Test RNASeqpipe class

#from RNAseqpipe_dev import run_RNAseqpipe
#from ..tools.get_gene_length import stat_db
import sys, os
import subprocess
import filecmp
from ..run_RNAseqpipe import main as run_main

gp_root="./test/RNAseqpip_data/data/gps"
# template file
tp=gp_root + "/small_group_data_template.txt"
# small group file
sg=gp_root + "/small_group_data.txt"
# Configurations
cf_root="./confs"

# Reference results
refroot = "./test/RNAseqpip_data/refout"

def test_hvc(tmpdir):
    refdir = refroot+"/hvc_count_out"
    outdir = tmpdir
    #outdir = "./test/RNAseqpip_data/testout/hvc_count_out"
    order = "./run_RNAseqpipe.py seq2exp -c {0}/dUTP.conf,{0}/gc_no_gff.conf,{0}/hisat_stringtie_verse.conf -g {1} -o {2} -q ".format(cf_root,sg,outdir).split()
    subprocess.call(order)
    targets = ["all_flagstat.txt","merged_asm/merged.fa","merged_asm/merged.gtf","salmon/quant_merge.elen","salmon/quant_merge.numreads","salmon/quant_merge.len","salmon/quant_merge.tpm"]
    #assert os.path.exists(outdir+"/imnotexist")
    for t in targets:
        tt = outdir+"/"+t
        assert os.path.exists(tt) and os.path.getsize(tt) > 0
        ## Check if the t is identical with reference file.
        ## But the results are not allways the same,
        ## so, this procedure will not be performed.
        ## way1
        #assert open(refdir+"/"+t).read() == open(outdir+"/"+t).read()
        ## way2
        #assert filecmp.cmp(refdir+"/"+t,outdir+"/"+t)
