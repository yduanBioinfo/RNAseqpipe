#!/usr/bin/env sh

input=$1

wd=`pwd`

data="${wd}/data"
# init conf
cat ${input} | awk -v data="${data}" '{gsub(/!datablock!/,data,$0);print}'
