#!/usr/bin/env bash

main_prog="../../../funcAnnot/b2gprog/blast2go.py"
genome="/biodb/genomes/Fish/grass_carp/new_genome/del_empty_line.alter.C_idella_female_scaffolds.fasta"
db1=/home/yduan/data/refseq/D_rerio/zebrafish.1.protein.faa
db2=/home/yduan/data/refseq/H_sapiens/human.protein.faa
${main_prog} --infile data/test.fa -o ./b2g_out -d ${db1} ${db2}
${main_prog} --infile data/test.gtf -s ${genome} -o ./b2g_gtf_out -d ${db1} ${db2}

# test_identity
for i in b2g_gtf_out/*
do
    diff ${i} test_out/${i##*/}
    if [[ $? -ne 0 ]];then
        echo ${i}
        exit 1
    fi
done

for i in b2g_out/*
do
    diff ${i} test_out/${i##*/}
    if [[ $? -ne 0 ]];then
        echo ${i}
        exit 1
    fi
done
rm -r b2g_out b2g_gtf_out
