# RNAseqpipe #

**A pipeline of RNAseq analysis**

Copyright (C) 2019 You Duan

## Prerequisites ##

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

本包是一个针对bulk-RNAseq数据分析的通用pipeline。整合了较为常用的Hisat-StringTie分析流程，可用于转录组表达量的定量工作和差异分析（该功能未完全实现）。  
本包建立了配置文件系统、元数据文件和多种分析流程。可以较为灵活地选择软件和配置参数。
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
- 数据系统(group data)
- 配置文件(Confs)
- 流程

### 一个简单的例子 ###

*输入测序数据，获得每个样本中每个基因的表达量*

```
# 建立RNAseqpipe_data文件夹。 
post_RNAseqpipe_install -o RNAseqpipe_data

# 切换目录
cd RNAseqpipe_data/test/test_bash/

# 运行示例程序
run_RNAseqpipe.py seq2exp -c data/confs/gc_no_gff.conf,../../confs/hisat_stringtie_verse.conf -g data/gps/small_group_data.txt -o testout/hvc_count_out
```

*各参数的含义*

-g 后面接[`数据文件`](#Group-data)  
gc_no_gff.conf为物种配置文件，定义了物种基因组文件的位置；  
hisat_stringtie_verse.conf为流程配置文件，定义了具体使用的软件和相应参数。  

-c 后面接[`配置文件`](#Confs)  
-o 后面接[`输出文件夹`](#Output)  
seq2exp为获取表达量子程序，相应的还有其他功能的子程序。  

```
# 通过帮助查看所有的子程序
$ run_RNAseqpipe.py --help
usage: run_RNAseqpipe.py [-h] {all,cptDE,seq2exp,func,verse,salmon}

RNA-seq analyse pipeline

positional arguments:
  {all,cptDE,seq2exp,func,verse,salmon}
                        all for whole pip/ali for alignment/ass for assembly/ cptDE for compute different expression

optional arguments:
  -h, --help            show this help message and exit
```

```
#通过帮助查看seq2exp子程序的全部参数
$ run_RNAseqpipe.py seq2exp --help

usage: run_RNAseqpipe.py [-h] [-1 [FQ1 [FQ1 ...]]] [-2 [FQ2 [FQ2 ...]]] [-c [CONF]] [-g [GROUP_DATA]] [-v] [-q] [-o [OUTPATH]]

RNA-seq seq2exp program

optional arguments:
  -h, --help            show this help message and exit
  -1 [FQ1 [FQ1 ...]], --fq1 [FQ1 [FQ1 ...]]
                        fq_1
  -2 [FQ2 [FQ2 ...]], --fq2 [FQ2 [FQ2 ...]]
                        fq_2
  -c [CONF], --conf [CONF]
                        configuration file
  -g [GROUP_DATA], --group_data [GROUP_DATA]
                        group_data file. conflict with -1 -2
  -v, --verbose         Out put all running information. Typically used in debug.
  -q, --quite           Running quitely.
  -o [OUTPATH], --outpath [OUTPATH]
                        outpath

```

### Group data ###
*数据文件系统*  
以配置文件的方式接收输入文件，这样做的好处首先是可以更好记录运行输入文件，其次可以自定义输出文件的名字，更为重要的是，其可以输入数据相关的表型值，可用于差异分析（虽然有设计，但未在实际的流程中应用）。
定义数据的配置文件可以在"data/gps"文件夹中找到。
以文件 **data/gps/small_group_data.txt** 为例

```
<dataset1>
filetype=fastq
library=PE#SE
TEST1   /home/yduan/scripts/python/projects/RNAseqpipe/test/RNAseqpip_data/data/fqs/test1_1.fq  /home/yduan/scripts/python/projects/RNAseqpipe/test/RNAseqpip_data/data/fqs/test1_2.fq
TEST2   /home/yduan/scripts/python/projects/RNAseqpipe/test/RNAseqpip_data/data/fqs/test2_1.fq  /home/yduan/scripts/python/projects/RNAseqpipe/test/RNAseqpip_data/data/fqs/test2_2.fq
TEST3   /home/yduan/scripts/python/projects/RNAseqpipe/test/RNAseqpip_data/data/fqs/test3_1.fq  /home/yduan/scripts/python/projects/RNAseqpipe/test/RNAseqpip_data/data/fqs/test3_2.fq
TEST4   /home/yduan/scripts/python/projects/RNAseqpipe/test/RNAseqpip_data/data/fqs/test3_1.fq  /home/yduan/scripts/python/projects/RNAseqpipe/test/RNAseqpip_data/data/fqs/test3_2.fq
</dataset1>

<dataset2>
filetype=fastq
library=SE
</dataset2>

<GROUP>
sample  bodysize
TEST1   Big
TEST2   Big
TEST3   Small
TEST4   Small
</GROUP>
```

group_data.txt使用xml格式，分为两个部分，数据（dataset）和分组（GROUP）。
**dataset**标签内定义了测序数据的基本属性和位置信息。
每一行从左至右为：数据的标签（样本名字），reads1文件位置，[reads2文件位置]（单端数据没有reads2）

**GROUP**标签定义了数据的分组信息，甚至可以包含多列的表型信息，是差异分析中需要使用到的信息。

### Confs ###
*配置文件系统*  
配置文件里面包含运行软件的路径（不提供具体路径时，使用环境变量中的软件），具体参数，以及部分必须文件的文件路径（例如基因组文件以及基因组注释文件gff）。
配置文件放在**confs**文件夹中。
以data/confs/gc_no_gff.conf（用于草鱼的转录组分析配置文件）为例：

```
<all>
genome=/home/yduan/scripts/python/projects/RNAseqpipe/test/RNAseqpip_data/data/test_genome.fa
</all>

<hisat2>
-x=/home/yduan/scripts/python/projects/RNAseqpipe/test/RNAseqpip_data/data/test_genome_idx
</hisat2>
```
    
*<all>* 标签中定义了基因组文件位置，还可定义注释数据库位置，总线程数等参数。
*<hisat2>* 标签中定义的是hisat2软件运行中使用的参数。
#### 参数优先级 ####

在简单的例子中，我们使用了2个配置文件，但配置文件系统远比这个灵活。  
我们在下面示例命令中，使用了三个配置文件（conf0,conf1,conf）：

    ${RNApipe} seq2exp -c ${conf0}","${conf1}","${conf} -g ${group_data} -o ${outdir}
    
理论上，RNAseqpipe可以接收任意多的配置文件，在本例中，其中conf的优先级最高，其次conf1，再其次conf0。意味着，如果三个文件中同时定义了同一参数，那么最终使用的参数为优先级最高的文件中定义的参数。
多参数系统事实上是为了适应多变的RNAseq测序手段和较为固定的策略的妥协方案。不同的转录组测序对应了不同的软件运行参数，因此参数往往是多变的；而往往有些转录组测序方案较为固定，如dUTP测序，同一实验室往往只研究少数的几个物种，例如我们实验室更多研究草鱼，在这些策略下，参数又是相对固定的。  
将conf文件具体化后，我们可将代码改写：
```
run_RNAseqpipe.py seq2exp -c ${cf_root}/dUTP.conf,${conf_gc_no_gff},${cf_root}/hisat_stringtie_verse.conf -g ${group_data} -o testout/hvc_count_out
```
在以上示例中，conf0对应了dUTP测序参数，conf1对应了草鱼的转录组分析参数，conf对应了使用hisat_stringtie_verse 的分析策略的参数。

### 流程 ###

一个子程序可以包括多种流程。如seq2exp程序中，可以选择用StringTie作为组装软件（hisat_stringtie_verse），也可以采用cufflinks作为组装软件（hisat_cuff_verse）。  
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

## Output ##
*结果文件*

### 目录结构 ###
```
├── all_flagstat.txt    （比对情况统计）
├── merged22_02_07_14.count （比对计数）
├── merged_asm  （组装结果文件夹）
│   ├── merged.fa  （fasta文件，组装结果或者是输入文件，取决于具体流程） 
│   ├── merged.gtf  （gtf文件，组装结果或者是输入文件，取决于具体流程）
│   └── salmon_index    （基于merged.fa的salmon索引文件）
├── salmon  （salmon表达量定量结果）
│   ├── quant_merge.elen （每个转录本有效长度）
│   ├── quant_merge.len （每个转录本长度）
│   ├── quant_merge.numreads    （每个转录本的表达量计数）
│   ├── quant_merge.tpm （每个转录本的TPM表达量）
│   ├── TEST1   （salmon计数的中间文件）
│   ├── TEST2   （salmon计数的中间文件）   
│   ├── TEST3   （salmon计数的中间文件）
│   └── TEST4   （salmon计数的中间文件）
├── TEST1   (样本TEST1的中间文件)
│   ├── flagstat.txt    （比对统计信息）
│   ├── sort.bam    （比对文件，已经排序）
│   ├── sort.exon.summary.txt   （计数信息）
│   ├── sort.exon.txt   （计数信息）
│   └── transcripts.gtf （组装文件）
├── TEST2   (样本TEST2的中间文件)
├── TEST3   (样本TEST3的中间文件)
├── TEST4   (样本TEST4的中间文件)
```
### File format ###
***all_flagstat.txt***

该文件是依赖`samtools flagstat`命令创建的，详细说明可参考[对应文档](http://www.htslib.org/doc/samtools-flagstat.html)

## to-do ##

repair logging from subprocess.communicate().

**Contact**

yduan94@ihb.ac.cn  
duanyou@outlook.com  