#!/bin/bash

if [ $# -ne 1 ];
then
    echo "Usage: 
        job_input.sh test.cfg"
    exit 1
fi

if [ ! -f $1 ]
then 
  echo "config file does not exist"
else 
  source $1

mkdir ${prefix} && cd ${prefix}
mkdir reads trimmedReads assembledGenome standardizedGenome annotatedGenome

WORKDIR=$(pwd)
read1_novo=${WORKDIR}/trimmedReads/forward_readP.fastq
read2_novo=${WORKDIR}/trimmedReads/reverse_readP.fastq

cat ${path_to_repo}/config_novo.txt |sed "s|WORKDIR|$WORKDIR|g"| sed "s|test|$prefix|" |sed "s|max_memory|${max_memory}|" | \
sed "s|path_to_seed|${seed}|" |sed "s|path_to_reference|${ref_fasta}|" |\
sed "s|read1|$read1_novo|" | sed "s|read2|$read2_novo|" > ${prefix}_config.txt

if [ $(echo "${read1##*.}") == "gz" ]
then 
cat ${path_to_repo}/Snakefile_single.py > _tmp.${prefix}_snakefile.py
else 
cat ${path_to_repo}/Snakefile_single.py |sed "/rule decompress/,/rule trimmomatic/{//p;d;}"|grep -v "rule decompress"| \
sed "s|WORKDIR/reads/forward_read.fastq|read1|g" |sed "s|WORKDIR/reads/reverse_read.fastq|read2|g" > _tmp.${prefix}_snakefile.py
fi

cat _tmp.${prefix}_snakefile.py |sed "s|WORKDIR|$WORKDIR|g"|sed "s|read1|$read1|g" | \
sed "s|read2|$read2|g" |sed "s|threads_available|$threads_available|" |sed "s|prefix|$prefix|g" |\
sed "s|path_to_ref_fasta|${ref_fasta}|g" |sed "s|path_to_ref_gb|${ref_gb}|g" |\
sed "s|path_to_trimmomatic|${path_to_trimmomatic}|g" |sed "s|path_to_novoplasty|${path_to_novoplasty}|g" |\
sed "s|path_to_repo|${path_to_repo}|g"|sed "s|path_to_PGA|${path_to_PGA}|g"  > ${prefix}_snakefile.py

rm _tmp.${prefix}_snakefile.py

#running snakemake 
snakemake --unlock
snakemake -s ${prefix}_snakefile.py -j ${threads_available}
fi
