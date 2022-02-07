#!/usr/bin/env sh

input=$1

wd=`pwd`

# init fqs
data_root="${wd}/data/fqs"
data="TEST1	${data_root}/test1_1.fq	${data_root}/test1_2.fq
TEST2	${data_root}/test2_1.fq	${data_root}/test2_2.fq
TEST3	${data_root}/test3_1.fq	${data_root}/test3_2.fq
TEST4	${data_root}/test3_1.fq	${data_root}/test3_2.fq"

# init gps
cat ${input} | awk -v data="${data}" '{if($0~/!datablock!/){print(data)}else{print}}'
