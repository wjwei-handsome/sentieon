[TOC]



# Sentieon SNP Calling

## 1.介绍

- ### 如果对下面的参数有什么不懂的地方，可以参考Sentieon官方文档

## 2.运行示例

### 	1.对所有样本进行单独的SNP Calling

```bash
bash DNAscope_gvcf.sh -s your_sample_name -q fastq_file1 -Q fastq_file2 -r your_ref_file -n threads -k(option) yes
```

​	**参数说明**：

| 参数名称 |            参数详解             |
| :------: | :-----------------------------: |
|    -s    |          你的样本名称           |
|    -q    |          双端测序之一           |
|    -Q    |          双端测序之二           |
|    -r    |         参考基因组文件          |
|    -n    |            线程数量             |
|    -k    | 是否保留测序文件，保留则键入yes |

### 	2.结果说明

​	**会为每个样本产生单独的结果文件夹，其中包含日志文件、gvcf结果已经各种统计信息**

### 	3.对所有样品的gvcf进行joint calling

```bash
bash joint_call_gvcfs.sh -r your_ref_file -o output_file
```

​	**参数说明**：

| 参数名称 |    参数详解    |
| :------: | :------------: |
|    -r    | 参考基因组文件 |
|    -o    |  输出结果文件  |

### 	4.结果说明

​	**会根据你的设定产生相应的结果文件**

## 3.脚本内容

**1.DNAscope_gvcf.sh**

```shell 
#!/usr/bin/env bash

USAGE="
---------------------------------------------------------------------------------------------------------------------                                                  
Sentieon  ---------  extremely fast SNP Calling shell script   
---------------------------------------------------------------------------------------------------------------------                                           
USAGE:
    bash $(basename $0) -s your_sample_name -q fastq_file1 -Q fastq_file2 -r your_ref_file -n threads -k(option) yes
---------------------------------------------------------------------------------------------------------------------                                  
"
### arguments setting ###
while getopts ":s:q:Q:r:n:k:" opt
do
    case $opt in
        s)
        sample_name=$OPTARG
        ;;
        q)
        RD1=$OPTARG
        ;;
        Q)
        RD2=$OPTARG
        ;;
        r)
        REF=$OPTARG
        ;;
        n)
        threads=$OPTARG
        ;;
        k)
        IF_rm_fastq=$OPTARG
        ;;
        ?)
        echo ${USAGE}
        exit 1;;
    esac
done


### check if file exist ###
if [ ! -f $RD1 ];then
    echo "Error:Fastq1 File ${RD1} not exist ";exit 
fi

if [ ! -f $RD2 ];then
    echo "Error:Fastq2 File ${RD2} not exist ";exit
fi

if [ ! -f $REF ];then
    echo "Error:Reference File $REF not exist ";exit
fi

### create work folder per sample ###
workdir=`pwd`/${sample_name}
[ ! -d $workdir ] && mkdir -p $workdir
cd $workdir

### load license ### ***PLEASE repalce the soft and license in your cluster***
export SENTIEON_LICENSE=mn01:9000
SENTIEON_INSTALL_DIR=/public/home/software/opt/bio/software/Sentieon/201808.07
### read group info for BWA ### (Options for platform ILLUMINA)
RG_info="@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA"

### standard out in logfile ####
logfile=${sample_name}_run.log
exec >$logfile 2>&1

####Sentieon proprietary compression
bam_option="--bam_compression 1"

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 

# speed up memory allocation malloc in bwa
export LD_sample_nameLOAD=$SENTIEON_INSTALL_DIR/lib/libjemalloc.so
export MALLOC_CONF=lg_dirty_mult:-1

( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R ${RG_info} -t $threads -K 10000000 $REF $RD1 $RD2 || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort $bam_option -r $REF -o ${sample_name}_sorted.bam -t $threads --sam2bam -i -
# options for move or kepp
if [ $IF_rm_fastq ];then
    echo "keep the fastq files"
else:
    rm ${RD1};rm ${RD2}
fi 
# ******************************************
# 2. Metrics
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $REF -t $threads -i ${sample_name}_sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o gc-report.pdf gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qd-report.pdf qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o mq-report.pdf mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o is-report.pdf is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads (Options)
# To mark duplicate reads only without removing them, remove "--rmdup" in the second command
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $threads -i ${sample_name}_sorted.bam --algo LocusCollector --rmdup --fun score_info score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $threads -i ${sample_name}_sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt $bam_option ${sample_name}_deduped.bam

# ******************************************
rm ${sample_name}_sorted.bam
# ******************************************
# 5. Base recalibration  BQSR
# ******************************************

# Perform recalibration
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $REF -t $threads -i ${sample_name}_deduped.bam --algo QualCal recal_data.table

# Perform post-calibration check (optional)
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $REF -t $threads -i ${sample_name}_deduped.bam -q recal_data.table --algo QualCal recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $threads --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv   
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o recal_plots.pdf recal.csv

# DNAscope
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $REF -t $threads -i ${sample_name}_deduped.bam -q recal_data.table --algo DNAscope  --emit_mode gvcf --emit_conf 10 --call_conf 10  ${sample_name}_scope.gvcf.gz
###这一步可以不做 但是做了万无一失
bgzip -d -@ $threads ${sample_name}-scope.gvcf.gz
rm ${sample_name}-scope.gvcf.gz.tbi
$SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert ${sample_name}_scope.gvcf ${sample_name}_scope.gvcf.gz
```

**2.joint_call_gvcfs.sh**

```shell
#!/usr/bin/env bash

##load license
export SENTIEON_LICENSE=mn01:9000
SENTIEON_INSTALL_DIR=/public/home/software/opt/bio/software/Sentieon/201808.07

UASGE="USAGE:  bash $(basename $0) -r your_ref_file -o output_file "

## arguments setting
while getopts ":r:o:" opt
do
    case $opt in
        r)
        REF=$OPTARG
        ;;
        o)
        output_file=$OPTARG
        ;;
        ?)
        echo ${USAGE}
        exit 1;;
    esac
done

##check arguments and file
if [ ! $REF ];then
    echo $UASGE;exit 
fi
if [ ! $output_file ];then
    echo $UASGE;exit 
fi
if [ ! -f $REF ];then
    echo "Error:Reference File $REF not exist ";exit
fi
## create a gvcf files list
find ./ -name '*_scope.gvcf.gz' > gvcfs.list
## 循环写变量
gvcf_argument=''
while read LINE;do
    gvcf_argument=${gvcf_argument}" -v ${LINE}"
done < gvcfs.list
## remove file list
rm gvcfs.list
## joint call all gvcfs
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r ${REF} --algo GVCFtyper ${gvcf_argument} ${output_file}
```

## 4.备注说明

1.sentieon中除了DNAscope，还有DNA-seq另一种算法，具体差异可以参考https://www.biorxiv.org/lookup/doi/10.1101/396325（Computational performance and accuracy of Sentieon DNASeq variant calling workflow)

2.其他具体参数可以参考官方文档