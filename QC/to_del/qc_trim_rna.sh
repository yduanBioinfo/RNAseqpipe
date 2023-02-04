#!/usr/bin/env bash

version_info="0.01"
help_info="
`basename "$0"`

This program is used to control quality of RNA-seq data.

License: GNU General Public License v2.0 (http://www.gnu.org/licenses/gpl-3.0.html)
Author: Mr. You Duan
Email: duanyou@outlook.com

Usage:
    `basename "$0"`  [options]  parameters

Options:
    -h, --help:     Display this message, then exit
    -v, --version:  Print version information
    -i:     Input file. If not provided, pop the first value in
        {parameters}, or STDIN
    -o:     Output file. If not provided, pop the first value in
        {parameters}, or STDOUT

parameters:
    Input file, Output file in order, if some of these augments are not provided as options.

"
#set -o noexec
#set -o verbose
#set -o xtrace

script_dir=$(dirname $(realpath $0))

exitInfo() {
    code=$1
    shift
    echo -e $@
    exit ${code}
}

if [[ -z "$1" ]]; then
    echo -e "${help_info}"
    exitInfo 0
fi

# set up option format
cmdnm=`basename "$0"`
# getopt -o hvi:o: -l help -l version: -- --help 4433 --version=3 asfsad
# --help --version '3' -- '4433' 'asfsad'
args=`getopt -o hvi:o: -l cores: -l -l help -l version -- "$@"`
eval set -- "$args"
# get options
for i in `seq $#`
do
    if [ $i -gt $# ]; then
        break
    fi
    case "${!i}" in
        -h |--help) echo -e "${help_info}"
            exitInfo 0;;
        -v |--version) exitInfo 0 "${version_info}";;
        -i) shift; fsrc="${!i}";;
        -o) shift; fobj="${!i}";;
        --) ;;
        *) if [[ -z "${fsrc}" ]]; then 
                fsrc="${!i}"; 
           elif [[ -z "${fobj}" ]]; then
                fobj="${!i}"; 
           else 
               echo -e "${help_info}" 
               exitInfo 1 "There are errors in options!"
           fi ;;
    esac
done

cores=1
outdir=filt_
clndir=clean
IlluQC_PRLL=/home/yduan/soft/bioinformatics/omics/QC/NGSQCToolkit-2.3/QC/IlluQC_PRLL.pl
AmbiguityFiltering=/home/yduan/soft/bioinformatics/omics/QC/NGSQCToolkit-2.3/Trimming/AmbiguityFiltering.pl
adaptor=/home/yduan/soft/bioinformatics/omics/QC/NGSQCToolkit-2.3/QC/adaptor_lib/v1_v2_LT

if [ ! -d ${clndir} ];then
    mkdir ${clndir}
fi

if [ ! -d ${outdir} ];then
    mkdir ${outdir}
fi

while [[ $1 ]];do 
    ${IlluQC_PRLL} -pe $1 $2 $adaptor A -c $cores -o ${outdir}
    base1=`basename $1`
    base2=`basename $2`
    qc1=${outdir}/${base1}_filtered
    qc2=${outdir}/${base2}_filtered
    ${AmbiguityFiltering} -i $qc1 -irev $qc2 -t3 -t5
    rm ${qc1} ${qc2}
    mv ${qc1}_trimmed ${clndir}/${base1%%.*}.clean_read1.fq
    mv ${qc2}_trimmed ${clndir}/${base2%%.*}.clean_read2.fq
    gzip ${clndir}/${base1%%.*}.clean_read1.fq
    gzip ${clndir}/${base2%%.*}.clean_read2.fq
    shift
    shift
done
