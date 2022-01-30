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

mkdir testout
## set2exp count mode
## quite mode
run_RNAseqpipe.py seq2exp -c ${cf_root}/dUTP.conf,${cf_root}/gc_no_gff.conf,${cf_root}/hisat_stringtie_verse.conf -g ${sg} -o testout/hvc_count_out -q && echo "hsv seq2exp done!" || echo "hsv seq2exp faild!"
## verbose mode
#run_RNAseqpipe.py seq2exp -c ${cf_root}/dUTP.conf,${cf_root}/gc_no_gff.conf,${cf_root}/dUTP.conf,${cf_root}/hisat_stringtie_verse.conf -g ${sg} -o testout/hvc_count_out --verbose && echo "hsv seq2exp done!" || echo "hsv seq2exp faild!"
#echo "../../run_RNAseqpipe.py seq2exp -c ${cf_root}/hisat_stringtie_verse.conf -g ${sg} -o testout/hvc_count_out --verbose" 

## hisat cufflinks verse
## quite mode
run_RNAseqpipe.py seq2exp -c ${cf_root}/gc_no_gff.conf,${cf_root}/hisat_cuff_verse.conf -g ${sg} -o testout/hsv_out -q && echo "hcv seq2exp done!" || echo "hcv seq2exp faild!"
# verbose mode
#run_RNAseqpipe.py seq2exp -c ${cf_root}/gc_no_gff.conf,${cf_root}/hisat_cuff_verse.conf -g ${sg} -o testout/hsv_out --verbose && echo "hcv seq2exp done!" || echo "hcv seq2exp faild!"

## hisat verse
## quite mode
run_RNAseqpipe.py seq2exp -c ${cf_root}/gc.conf,${cf_root}/hisat_verse_count.conf -g ${sg} -o testout/hv_out -q && echo "hv seq2exp done!" || echo "hcv seq2exp faild!"
## verbose mode
#run_RNAseqpipe.py seq2exp -c ${cf_root}/gc.conf,${cf_root}/hisat_verse_count.conf -g ${sg} -o testout/hv_out --verbose && echo "hv seq2exp done!" || echo "hcv seq2exp faild!"

## salmon
run_RNAseqpipe.py salmon -g $sg --gtf data/dlmrna.gc.final.gff --genome /biodb/genomes/Fish/grass_carp/new_genome/del_empty_line.alter.C_idella_female_scaffolds.fasta -o testout/salmon && echo "salmon done!" || echo "salmon faild!"
