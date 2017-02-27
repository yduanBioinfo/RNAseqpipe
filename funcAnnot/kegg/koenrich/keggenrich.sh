#!/usr/bin/env sh

#[Usage]: ./keggenrich.sh kegg_annot_file

#./keggCountal.py $1 kodb/oneGcGene_oneKO kodb/pathway.ko $1".kg"
./keggCountal.py $1 kodb/keggannot_15_12_30.txt kodb/pathway.ko $1".kg"
./mykeggenrich.R $1".kg"
