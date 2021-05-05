# Test RNASeqpipe class

#from RNAseqpipe_dev import run_RNAseqpipe
#from ..tools.get_gene_length import stat_db
import sys, os
import subprocess
import filecmp
import glob
from ..run_RNAseqpipe import main as run_main

gp_root="./test/RNAseqpip_data/data/gps"
# template file
tp=gp_root + "/small_group_data_template.txt"
# small group file
sg=gp_root + "/small_group_data.txt"
# Configurations
cf_root="./confs"

# Requeired data
genome = "/biodb/genomes/Fish/grass_carp/new_genome/del_empty_line.alter.C_idella_female_scaffolds.fasta"

# Reference results
refroot = "./test/RNAseqpip_data/refout"

def _check_exist(targets):
    for tt in targets:
        print("-------------------------Here is tt-----------------------")
        print(tt)
        assert os.path.exists(tt) and os.path.getsize(tt) > 0

def test_hvc(tmpdir):
    # Hisat2 + stringtie + verse + salmon
    #refdir = refroot+"/hvc_count_out"
    order = "./run_RNAseqpipe.py seq2exp -c {0}/dUTP.conf,{0}/gc_no_gff.conf,{0}/hisat_stringtie_verse.conf -g {1} -o {2} ".format(cf_root,sg,tmpdir).split()
    subprocess.call(order)
    
    targets = ["all_flagstat.txt","merged_asm/merged.fa","merged_asm/merged.gtf","salmon/quant_merge.elen","salmon/quant_merge.numreads","salmon/quant_merge.len","salmon/quant_merge.tpm"]
    for t in targets:
        tt = tmpdir.join(t)
        assert os.path.exists(tt) and os.path.getsize(tt) > 0
    #    ## Check if the t is identical with reference file.
    #    ## But the results are not allways the same,
    #    ## so, this procedure will not be performed.
    #    ## way1
    #    #assert open(refdir+"/"+t).read() == open(tmpdir.join(t)).read()
    #    ## way2
    #    #assert filecmp.cmp(refdir+"/"+t,tmpdir.join(t))

def test_hcv(tmpdir):
    # Hisat2 + cufflinks + verse
    #refdir = refroot+"/hcv_out"
    order = "./run_RNAseqpipe.py seq2exp -c {0}/gc_no_gff.conf,{0}/hisat_cuff_verse.conf -g {1} -o {2} ".format(cf_root,sg,tmpdir).split()
    subprocess.call(order)
    
    targets = ["all_flagstat.txt","merged_asm/merged.fa","merged_asm/merged.gtf","salmon/quant_merge.elen","salmon/quant_merge.numreads","salmon/quant_merge.len","salmon/quant_merge.tpm"]
    for t in targets:
        tt = tmpdir.join(t)
        assert os.path.exists(tt) and os.path.getsize(tt) > 0

def test_hv(tmpdir):
    # Hisat2 + verse
    #refdir = refroot+"/hcv_out"
    order = "./run_RNAseqpipe.py seq2exp -c {0}/gc.conf,{0}/hisat_verse_count.conf -g {1} -o {2}".format(cf_root,sg,tmpdir).split()
    subprocess.call(order)
    
    targets = ["all_flagstat.txt"]
    targets = map(tmpdir.join, targets)
    _check_exist(targets)
    #for t in targets:
    #    tt = tmpdir.join(t)
    #    assert os.path.exists(tt) and os.path.getsize(tt) > 0
    # for merged*.count

    _check_exist(glob.glob(str(tmpdir.join("merged*.count"))))


def test_salmon(tmpdir):
    # salmon
    #refdir = refroot+"/hcv_out"
    gff = "./test/RNAseqpip_data/data/dlmrna.gc.final.gff"
    order = "./run_RNAseqpipe.py salmon --gtf {0} -g {1} -o {2} --genome {3}".format(gff,sg,tmpdir,genome).split()
    subprocess.call(order)
    
    targets = ["all_flagstat.txt", "quant_merge.elen","quant_merge.numreads","quant_merge.len","quant_merge.tpm"]
    for t in targets:
        tt = tmpdir.join(t)
        assert os.path.exists(tt) and os.path.getsize(tt) > 0

## hisat cufflinks verse
## quite mode
#../../run_RNAseqpipe.py seq2exp -c ${cf_root}/gc_no_gff.conf,${cf_root}/hisat_cuff_verse.conf -g ${sg} -o testout/hsv_out -q && echo "hcv seq2exp done!" || echo "hcv seq2exp faild!"
# verbose mode
#../../run_RNAseqpipe.py seq2exp -c ${cf_root}/gc_no_gff.conf,${cf_root}/hisat_cuff_verse.conf -g ${sg} -o testout/hsv_out --verbose && echo "hcv seq2exp done!" || echo "hcv seq2exp faild!"

## hisat verse
## quite mode
#../../run_RNAseqpipe.py seq2exp -c ${cf_root}/gc.conf,${cf_root}/hisat_verse_count.conf -g ${sg} -o testout/hv_out -q && echo "hv seq2exp done!" || echo "hcv seq2exp faild!"
## verbose mode
#../../run_RNAseqpipe.py seq2exp -c ${cf_root}/gc.conf,${cf_root}/hisat_verse_count.conf -g ${sg} -o testout/hv_out --verbose && echo "hv seq2exp done!" || echo "hcv seq2exp faild!"

## salmon
#../../run_RNAseqpipe.py salmon -g $sg --gtf data/dlmrna.gc.final.gff --genome /biodb/genomes/Fish/grass_carp/new_genome/del_empty_line.alter.C_idella_female_scaffolds.fasta -o testout/salmon

