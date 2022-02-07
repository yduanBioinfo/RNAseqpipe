# RNAseqpipe #

**A pipeline of RNAseq analysis**

Copyright (C) 2019 You Duan

general RNAseq process pipeline

**Contact**

yduan94@ihb.ac.cn

**Prerequisites**

Python > 3.6.7
hisat2 > 2.0.4
StringTie > 1.3.1c
verse > 0.1.5
htseq-count 0.6.1p1
salmon 
samtools

### python and R packages ###

#### R 包 ####
argparser
ggplot2
```
R
> install.packages('argparser')
> install.packages('ggplot2')
```

#### python 包 ####
```
# 可跳过 pyper
pip install pyper
pip install dypylib
```

**Introduction**

...

## Install ##
```
# 使用pip进行安装
pip install RNAseqpipe
# 在RNAseqpipe_data目录下生成配置文件和测试文件（RNAseqpipe_data
# 可换成任意文件名）。
post_RNAseqpipe_install -o RNAseqpipe_data
# 运行测试
cd RNAseqpipe_data/test/test_bash/
sh cmd
```
or
```
# 使用超级用户权限进行配置
# 使用pip进行安装
sudo pip install RNAseqpipe
# 在RNAseqpipe_data目录下生成配置文件和测试文件（RNAseqpipe_data
# 可换成任意文件名）。
sudo post_RNAseqpipe_install -o RNAseqpipe_data
# 运行测试
cd RNAseqpipe_data/test/test_bash/
sh cmd
```


## Usage ##

使用包括三方面内容：
- 数据系统
- 配置文件
- 流程

### 数据系统 ###
以配置文件的方式接收输入文件，这样做的好处首先是可以更好记录运行输入文件，其次可以自定义输出文件的名字，更为重要的是，其可以输入数据相关的表型值，可用于差异分析（虽然有设计，但未在实际的流程中应用）。
定义数据的配置文件可以在"group_data"文件夹中找到。
以文件 **group_data/group_data.txt** 为例

    <dataset1>
    filetype=fastq
    library=PE#SE
    15_3_BB1        /home/yduan/data/growth/ruibo/2015-3-BB1.read1_Clean.fastq      /home/yduan/data/growth/ruibo/2015-3-BB1.read2_Clean.fastq
    15_3_BB2        /home/yduan/data/growth/ruibo/2015-3-BB2.read1_Clean.fastq      /home/yduan/data/growth/ruibo/2015-3-BB2.read2_Clean.fastq
    15_3_SB2        /home/yduan/data/growth/ruibo/2015-3-SB2.read1_Clean.fastq      /home/yduan/data/growth/ruibo/2015-3-SB2.read2_Clean.fastq
    15_3_SB5        /home/yduan/data/growth/ruibo/2015-3-SB5.read1_Clean.fastq      /home/yduan/data/growth/ruibo/2015-3-SB5.read2_Clean.fastq
    </dataset1>
    
    15_3_SB1        /home/yduan/data/growth/ruibo/2015-3-SB1.read1_Clean.fastq      /home/yduan/data/growth/ruibo/2015-3-SB1.read2_Clean.fastq
    15_3_SB1        Small
    <dataset2>
    filetype=fastq
    library=SE
    </dataset2>
    
    <GROUP>
    sample  bodysize
    15_3_BB1        Big
    15_3_BB2        Big
    15_3_SB2        Small
    15_3_SB5        Small
    </GROUP>

group_data.txt使用xml格式，分为两个部分，数据（dataset）和分组（GROUP）。
**dataset**标签内定义了测序数据的基本属性和位置信息。
每一行位置信息可以包含：数据的标签，reads1位置，[reads2位置]（单端数据没有reads2）

**GROUP**标签定义了数据的分组信息，甚至可以包含多列的表型信息，是差异分析中需要使用到的信息。

### 配置文件系统 ###
配置文件里面包含运行软件的路径（不提供具体路径时，使用环境变量中的软件），具体参数，以及部分必须文件的文件路径（例如基因组文件以及基因组注释文件gff）。
配置文件放在**confs**文件夹中。
以confs/gc.conf（用于草鱼的转录组分析配置文件）为例：

    <all>
    #general parameters
    genome=/home/lab/genomes/Fish/grass_carp/new_genome/del_empty_line.alter.C_idella_female_scaffolds.fasta
    #gff=/home/lab/genomes/Fish/grass_carp/new_genome/dlmrna.gc.final.gff
    gff=/home/lab/genomes/Fish/grass_carp/new_genome/dlmrna.gc.final.gtf
    kodb=/home/yduan/data/growth/clean/polya/analysis/merged_asm/koannot/merged.koannot#ke gg annot db
    godb=/home/yduan/data/growth/clean/polya/analysis/merged_asm/goannot/merged_all.annot
    kgmap=/home/yduan/data/funcdb/kodb/pathway.ko#keggmap
    p=8
    </all>
    
    <hisat2>
    -x=/home/lab/genomes/Fish/grass_carp/new_genome/daC_ide
    </hisat2>
    
*<all>* 标签中定义了包括文件位置，注释数据库位置，总线程数等参数。
*<hisat2>* 标签中定义的是hisat2软件运行中使用的参数。
#### 参数优先级 ####


在示例文件中，使用了三个运行参数（conf0,conf1,conf）

    ${RNApipe} seq2exp -c ${conf0}","${conf1}","${conf} -g ${group_data} -o ${outdir}
    
理论上，RNAseqpipe可以接收任意多的配置文件，在本例中，其中conf的优先级最高，其次conf1，再其次conf0。意味着，如果三个文件中同时定义了同一参数，那么最终使用的参数为优先级最高的文件中定义的参数。
多参数系统事实上是为了适应多变的RNAseq测序手段和较为固定的策略的妥协方案。不同的转录组测序对应了不同的软件运行参数，因此参数往往是多变的；而往往有些转录组测序方案较为固定，如dUTP测序，同一实验室往往只研究少数的几个物种，例如我们实验室更多研究草鱼，在这些策略下，参数又是相对固定的。
在我们的示例文件中，conf0对应了dUTP测序参数，conf1对应了草鱼的转录组分析参数，conf对应了使用hisat_stringtie_verse 的分析策略的参数。

### 流程 ###

在hisat_stringtie_verse.conf 文件中<all> 标签下pipe参数的值为hsv，代表了使用hisat2 + stringtie + verse 的分析流程。这样的预定义流程包括如下：
- hisat_stringtie_verse
- hisat_cuff_verse
- hisat_htseq_count （不组装转录本，直接计数）
- hisat_verse_count （不组装转录本，直接计数）


**hisat_stringtie_verse.conf** 文件

    cat hisat_stringtie_verse.conf
    
    <all>
    # hisat stringtie htseq-count pipe.
    # Strand specific library, dUTP methods.
    pipe=hsv
    </all>


## to-do ##


repair logging from subprocess.communicate().

