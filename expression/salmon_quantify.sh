#!/usr/bin/env bash

Usage="salmon_quantify.sh ingff genome outdir threads library1_name fq1_1 fq1_2 library2_name fq2_1 fq2_2 ..."
inputgff=$1
genome=$2
outdir=$3
threads=$4
mkdir ${outdir}
shift;shift;shift;shift

# Creat transcriptome and build index
transcriptome=${inputgff%/*}/merged.fa
salmon_index=${transcriptome%/*}/salmon_index
if [[ -d ${salmon_index} ]]
then
    echo "index exist"
else
    gffread ${inputgff} -g ${genome} -w ${transcriptome}
    salmon index -t ${transcriptome} -i ${salmon_index}
fi

declare -a names
index=1
## non-alignment-based
while [[ $1 ]];do
    sample_name=$1
    fq1=$2
    fq2=$3
    names[${index}]=${outdir}/${sample_name}
    salmon quant -p ${threads} -i ${salmon_index} -l A -1 ${fq1} -2 ${fq2} -o ${outdir}/${sample_name} --validateMappings
    ((index++))
    shift;shift;shift
done

## merge results
salmon quantmerge --column len --quants ${names[@]} -o ${outdir}/quant_merge.len
salmon quantmerge --column elen --quants ${names[@]} -o ${outdir}/quant_merge.elen
salmon quantmerge --column tpm --quants ${names[@]} -o ${outdir}/quant_merge.tpm
salmon quantmerge --column numreads --quants ${names[@]} -o ${outdir}/quant_merge.numreads

