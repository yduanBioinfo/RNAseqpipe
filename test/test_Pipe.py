# Test RNASeqpipe class

#from RNAseqpipe_dev import run_RNAseqpipe
#from ..tools.get_gene_length import stat_db
import sys
import subprocess
from ..run_RNAseqpipe import main as run_main

gp_root="./test/RNAseqpip_data/data/gps"
# template file
tp=gp_root + "/small_group_data_template.txt"
# small group file
sg=gp_root + "/small_group_data.txt"
# Configurations
cf_root="./confs"

def test_len():
    assert 1==1

def test_which():
    outdir = "./test/RNAseqpip_data/testout/hvc_count_out"
    order = "./run_RNAseqpipe.py seq2exp -c {0}/dUTP.conf,{0}/gc_no_gff.conf,{0}/hisat_stringtie_verse.conf -g {1} -o {2} -q ".format(cf_root,sg,outdir).split()
    subprocess.call(order)
    assert order == "mm"
    #&& echo "hsv seq2exp done!" || echo "hsv seq2exp faild!"
