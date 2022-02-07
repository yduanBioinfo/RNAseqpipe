#!/usr/bin/env sh

# INIT of group_data
# Directory for group data files.
gp_root="data/gps"
# template file
tp="${gp_root}/small_group_data_template.txt"
# small group file
sg="${gp_root}/small_group_data.txt"
./init_test_data.sh ${tp} > ${sg}

# Configuration files
cf_root="../../confs"
gc_confs="data/confs"
tconf_gc_no_gff="${gc_confs}/gc_no_gff_template.conf"
tconf_gc="${gc_confs}/gc_template.conf"
conf_gc_no_gff="${gc_confs}/gc_no_gff.conf"
conf_gc="${gc_confs}/gc.conf"
./init_test_conf.sh $tconf_gc_no_gff > $conf_gc_no_gff 
./init_test_conf.sh $tconf_gc > $conf_gc 

mkdir testout
## set2exp count mode
## quite mode
run_RNAseqpipe.py seq2exp -c ${cf_root}/dUTP.conf,${conf_gc_no_gff},${cf_root}/hisat_stringtie_verse.conf -g ${sg} -o testout/hvc_count_out -q && echo "hsv seq2exp done!" || echo "hsv seq2exp faild!"
## verbose mode
#run_RNAseqpipe.py seq2exp -c ${cf_root}/dUTP.conf,${conf_gc_no_gff},${cf_root}/dUTP.conf,${cf_root}/hisat_stringtie_verse.conf -g ${sg} -o testout/hvc_count_out --verbose && echo "hsv seq2exp done!" || echo "hsv seq2exp faild!"
#echo "../../run_RNAseqpipe.py seq2exp -c ${cf_root}/hisat_stringtie_verse.conf -g ${sg} -o testout/hvc_count_out --verbose" 

## hisat cufflinks verse
## quite mode
run_RNAseqpipe.py seq2exp -c ${conf_gc_no_gff},${cf_root}/hisat_cuff_verse.conf -g ${sg} -o testout/hsv_out -q && echo "hcv seq2exp done!" || echo "hcv seq2exp faild!"
# verbose mode
#run_RNAseqpipe.py seq2exp -c ${conf_gc_no_gff},${cf_root}/hisat_cuff_verse.conf -g ${sg} -o testout/hsv_out --verbose && echo "hcv seq2exp done!" || echo "hcv seq2exp faild!"

## hisat verse
## quite mode
run_RNAseqpipe.py seq2exp -c ${conf_gc},${cf_root}/hisat_verse_count.conf -g ${sg} -o testout/hv_out -q && echo "hv seq2exp done!" || echo "hcv seq2exp faild!"
## verbose mode
#run_RNAseqpipe.py seq2exp -c ${conf_gc},${cf_root}/hisat_verse_count.conf -g ${sg} -o testout/hv_out --verbose && echo "hv seq2exp done!" || echo "hcv seq2exp faild!"

## salmon
run_RNAseqpipe.py salmon -g $sg --gtf "data/test.dlmrna.gc.final.gff" --genome "data/test_genome.fa" -o testout/salmon && echo "salmon done!" || echo "salmon faild!"
