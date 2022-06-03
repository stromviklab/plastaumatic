#!/bin/bash

if [ $# -ne 1 ];then
  echo "Usage: 
      run_plastaumatic.sh test.cfg"
  exit 1
fi

if [ ! -f $1 ];then 
  echo "config file does not exist"
else 
  source $1

  for line in `cat $fof`;do 
    prefix=$(echo $line|cut -d ',' -f1)
    read1=$(echo $line|cut -d ',' -f2)
    read2=$(echo $line|cut -d ',' -f3)
    
    if [[ -f $read1 && -f $read2 ]];then 
      cd $out
      mkdir -p ${prefix} && cd ${prefix}

      WORKDIR=$(pwd)
      read1_novo=${WORKDIR}/trimmedReads/forward_readP.fastq
      read2_novo=${WORKDIR}/trimmedReads/reverse_readP.fastq

      cat ${path_to_repo}/config_novo.txt |sed "s|WORKDIR|$WORKDIR|g"| sed "s|test|$prefix|" |sed "s|max_memory|${max_memory}|" | \
      sed "s|path_to_seed|${seed}|" |sed "s|path_to_reference|${ref_fasta}|" |sed "s|range|$range|"|\
      sed "s|read1|$read1_novo|" | sed "s|read2|$read2_novo|" > ${prefix}_config.txt

      if [ $(echo "${read1##*.}") == "gz" ];then 
        cat ${path_to_repo}/Snakefile > _tmp.${prefix}_snakefile
      else 
        cat ${path_to_repo}/Snakefile |sed "/rule decompress/,/rule trimmomatic/{//p;d;}"|grep -v "rule decompress"| \
        sed "s|WORKDIR/reads/forward_read.fastq|read1|g" |sed "s|WORKDIR/reads/reverse_read.fastq|read2|g" > _tmp.${prefix}_snakefile
      fi

      cat _tmp.${prefix}_snakefile |sed "s|WORKDIR|$WORKDIR|g"|sed "s|read1|$read1|g" | \
      sed "s|read2|$read2|g" |sed "s|threads_available|$threads_available|" |sed "s|prefix|$prefix|g" |\
      sed "s|path_to_ref_fasta|${ref_fasta}|g" |sed "s|path_to_ref_gb|${ref_gb}|g" |\
      sed "s|path_to_trimmomatic|${path_to_trimmomatic}|g" |sed "s|path_to_novoplasty|${path_to_novoplasty}|g" |\
      sed "s|path_to_repo|${path_to_repo}|g"|sed "s|path_to_PGA|${path_to_PGA}|g"|sed "s|ADAPTERS|${adapters}|"  > snakefile

      rm _tmp.${prefix}_snakefile

      # modify the ISC script
      cat ${path_to_repo}/check_internal_stops.sh|sed "s|translateDna|${path_to_repo}/translateDna.pl|g" > ${prefix}_isc.sh 
      chmod u+x ${prefix}_isc.sh 

      #running snakemake 
      snakemake --unlock
      snakemake -j ${threads_available} -q
    else
      echo "The read file/files for the $prefix does not exist"
    fi 
  done
fi
