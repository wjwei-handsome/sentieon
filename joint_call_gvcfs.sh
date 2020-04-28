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