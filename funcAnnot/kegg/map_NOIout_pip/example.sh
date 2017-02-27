echo "ko00010" | ./map2ko.py ../kodb/pathway.ko | ./ko2gene.py ../kegg_12_28_annot.txt | ./gene2annot.py ../rpkm_all_15_12_30.txt > myout.txt
