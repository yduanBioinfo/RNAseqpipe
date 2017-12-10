#!/usr/bin/env sh

#convert soft-clipping to hard-clipping
#conv_clipping.sh infile outfile
#infile and out file should be sam file

infile=$1
outfile=$2
awk 'BEGIN {OFS="\t"} {if (substr($1,1,1)=="@") {print;next}; split($6,C,/[0-9]*/); split($6,L,/[SMDIN]/); if (C[2]=="S") {$10=substr($10,L[1]+1); $11=substr($11,L[1]+1)}; if (C[length(C)]=="S") {L1=length($10)-L[length(L)-1]; $10=substr($10,1,L1); $11=substr($11,1,L1); }; gsub(/[0-9]*S/,"",$6); print}' $1 > $2

