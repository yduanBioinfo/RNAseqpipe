#!/usr/bin/env sh

sd=`dirname $0`
sd=`realpath ${sd}`

Assembly_root=$1
#mappint_rate="${Assembly_root}/all_samples_flagstat.txt"
#merge_prog="${sd}/merge_flagstat.py"
## Get mapping rate for every sample with samtools.
#for i in ${Assembly_root}/*/sort.bam
#do
#    samtools flagstat ${i} > ${i%/*}/flagstat.txt
#done
## Merge mapping rate

${merge_prog} ${Assembly_root}/*/flagstat.txt > ${Assembly_root}/all_flagstat.txt
