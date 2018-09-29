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
# set2exp count mode
../../RNAseqpip.py seq2exp -c ${cf_root}/hisat_stringtie_verse.txt -g ${sg} -o testout/hvc_count_out -q && echo "hsv seq2exp done!" || echo "hsv seq2exp faild!"

# hisat cufflinks verse
../../RNAseqpip.py seq2exp -c ${cf_root}/configuration.txt -g ${sg} -o testout/hsv_out -q && echo "hcv seq2exp done!" || echo "hcv seq2exp faild!"
