#!/usr/bin/env bash

b_root="/home/wdye/mybiodb/software/ncbi-blast-2.8.0+/bin/"
wd=$(dirname `realpath $0`)
version_info="0.01"
help_info="
`basename "$0"`

This program is used to run blast2go annotation.

License: GNU General Public License v2.0 (http://www.gnu.org/licenses/gpl-3.0.html)
Author: Mr. You Duan
Email: duanyou12345@126.com

Usage:
    `basename "$0"`  [options]  parameters

Options:
    -h, --help:     Display this message, then exit
    -v, --version:  Print version information
    -f: --fasta:    If the format of input file is fasta, set this option.
        Xml is the default format.
    -i:     Input file. If not provided, pop the first value in
        {parameters}, or STDIN
    -o:     Output directory. If not provided, pop the first value in
        {parameters}, or STDOUT

parameters:
    Input file, Output file in order, if some of these augments are not provided as options.

"

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
args=`getopt -o hvfi:o: -l help -l version -l fasta -- "$@"`
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
        -f) fasta="fasta";;
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

xml2b2g()
{
    # Run b2g annotation start with blastout results (.xml).
    local input=$1
    local outdir=$2
    # b2g programe only annot xml with BLASTX 2.2.27+ flag.
    # A more flexible substitution is in need.
    cat ${input} | sed 's/BLASTX 2.8.0+/BLASTX 2.2.27+/' > ${outdir}/blast_out_27.xml
    # Split xml into small files in case of abortion for file size.
    ${wd}/split_blast_xml.py ${outdir}/blast_out_27.xml ${outdir}/pysplit
    for i in ${outdir}/pysplit*.xml;do
        b2g -in ${i} -annot -out ${i%.xml}b2gres
    done
    # merge
    cat ${outdir}/pysplit*b2gres.annot > ${outdir}/pysplit_all.annot
    rm ${outdir}/pysplit*.xml ${outdir}/pysplit*b2gres.annot
}

myb2g()
{
    # Run b2g annotation start with fasta.
    local input=$1
    local outdir=$2
    ${b_root}/blastx -query ${input} -db ${index} -outfmt 5 -evalue 1E-5 -max_target_seqs 200 -num_threads 128 -show_gis -out ${outdir}/blast_out.xml
    ###  further TO-DO.
    ### integret alter BLASTX in split_blast_xml.py header data.
    ### muilt-threading
    xml2b2g ${outdir}/blast_out.xml ${outdir}
    #cat ${outdir}/blast_out.xml | sed 's/BLASTX 2.8.0+/BLASTX 2.2.27+/' > ${outdir}/blast_out_27.xml
    #${wd}/split_blast_xml.py ${outdir}/blast_out_27.xml ${outdir}/pysplit
    #for i in ${outdir}/pysplit*.xml;do
    #    b2g -in ${i} -annot -out ${i%.xml}b2gres
    #done
    ## merge
    #cat ${outdir}/pysplit*b2gres.annot > ${outdir}/pysplit_all.annot
    #rm ${outdir}/pysplit*.xml ${outdir}/pysplit*b2gres.annot
}

annotation="${fobj}"
mkdir ${annotation}
blast_xml="${annotation}/merged2Tele.xml"
if [[ ${fasta} ]];then
    myb2g ${fsrc} ${annotation}
else
    xml2b2g ${fsrc} ${annotation}
fi


