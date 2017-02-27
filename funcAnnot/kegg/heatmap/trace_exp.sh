#plot heatmap of enriched komap

#[Usage]: ./heatmap.sh kegg.adjDiff expfile annotfile mapfile
prefix=`dirname $(readlink $0 || echo $0)`
#cat ../kegg_12_28_annot.txt.kg.adjDiff | awk 'BEGIN{expfile="../rpkm_all_15_12_30.txt";annotfile="../kegg_12_28_annot.txt";mapfile="~/growth/pip/KEGG/kodb/pathway.ko"}/^ko/{system("echo "$1 "|~/growth/pip/KEGG/map_NOIout_pip/map2ko.py "mapfile"| ~/growth/pip/KEGG/map_NOIout_pip/ko2gene.py "annotfile"| ~/growth/pip/KEGG/map_NOIout_pip/gene2annot.py "expfile "|./heatmap.R "$1)}'

#cat $1 | awk 'BEGIN{expfile='$2';annotfile='$3';mapfile='$4'}/^ko/{system("echo "$1 "|/home/yduan/script/python/RNAseqpip/funcAnnot/kegg/map_NOIout_pip/map2ko.py "mapfile"| /home/yduan/script/python/RNAseqpip/funcAnnot/kegg/map_NOIout_pip/ko2gene.py "annotfile"| /home/yduan/script/python/RNAseqpip/funcAnnot/kegg/map_NOIout_pip/gene2annot.py "expfile"|./heatmap.R "$1)}'

cat $1 | awk 'BEGIN{expfile="'$2'";annotfile="'$3'";mapfile="'$4'";prefix="'$prefix'"}/^ko/{system("echo "$1 "|/home/yduan/script/python/RNAseqpip/funcAnnot/kegg/map_NOIout_pip/map2ko.py "mapfile"| /home/yduan/script/python/RNAseqpip/funcAnnot/kegg/map_NOIout_pip/ko2gene.py "annotfile"| /home/yduan/script/python/RNAseqpip/funcAnnot/kegg/map_NOIout_pip/gene2annot.py "expfile)}'
