#!/bin/bash

## Exits if a command fails
set -eo pipefail

## parse arguments
while getopts ':s:g:r:o:f:n:h' options
do
  case $options in
    s) seed=$OPTARG ;;
    g) ref_gb=$OPTARG ;;
	r) range=$OPTARG ;;
	o) out=$OPTARG ;;
    f) fof=$OPTARG ;;
	n) novo=$OPTARG ;;
	h) echo -e "Usage:\tplastaumatic -s seed.fa -g reference.gb -r <120000-160000> -o out_dir -f fof.txt -n NOVOPlasty4.3.1.pl\
				\noptions:\n\t -s\t Path to the seed file for assembly\
				\n\t -g\t Path to the reference GenBank file\
				\n\t -r\t Plastome assembly size range [150000-160000]\
				\n\t -o\t Path to the output directory\
				\n\t -f\t Path to the file-of-filenames with reads\
				\n\t -n\t Path to NOVOPlasty executable\
				\n\t -h\t Shows this help message\n"
		exit;;
	:) echo "ERROR: -$OPTARG requires an argument" 1>&2
    	exit;;
	\?) echo "Invalid Option: try  plastaumatic -h" 1>&2
      exit;;
  esac
done

## Checks missing options
if [ $# -eq 0 ] || [ -z "$seed" ] || [ -z "$ref_gb" ] || [ -z "$range" ] || [ -z "$out" ] || [ -z "$fof" ]  || [ -z "$novo" ];then 
	echo "ERROR: missing options, try
            plastaumatic -h"
	exit
fi

## Checks if the input files are present
if [ ! -f $seed ];then
	echo -e "ERROR:\t Seed file does not exist, check again" && exit
elif [ ! -f $ref_gb ];then 
	echo -e "ERROR:\t Reference GenBank file does not exist, check again" && exit
elif ! command -v fastp &> /dev/null ;then 
	echo -e "ERROR:\t fastp executable does not exist in the path, check again" && exit
elif ! command -v samtools &> /dev/null ;then 
	echo -e "ERROR:\t samtools executable does not exist in the path, check again" && exit
elif ! command -v blastn &> /dev/null ;then 
	echo -e "ERROR:\t blastn executable does not exist in the path, check again" && exit
elif [ ! -f $novo ];then 
	echo -e "ERROR:\t NOVOPlasty executable does not exist, check again" && exit
elif [ ! -f $fof ];then 
	echo -e "ERROR:\t fof file does not exist, check again" && exit
fi

## Sets memory and cpus available
# threads=$(nproc --all)
# max_memory=$(awk '/MemAvailable/ {printf "%.0f\n", $2/1024/1024*0.9}' /proc/meminfo)
threads=1
max_memory=30
repo=$(dirname $0)
seed_path=$(pwd)/${seed#$(pwd)/}
ref_gb_path=$(pwd)/${ref_gb#$(pwd)/}
novo_path=$(pwd)/${novo#$(pwd)/}
out_path=$(pwd)/${out#$(pwd)/}

for line in `cat $fof`;do 
	prefix=$(echo $line|cut -d ',' -f1)
	read1=$(echo $line|cut -d ',' -f2)
	read2=$(echo $line|cut -d ',' -f3)
	base1=$(basename $read1)
	base2=$(basename $read2)

	if [ ! -f $read1 ] || [ ! -f $read2 ];then
		echo -e "ERROR: \t Read files does not exist, check again" && exit 1
	fi  

	echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tRunning plastaumatic on ${prefix}" 

	## Setting the directories
	echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tsetting up the directories" 
	mkdir -p ${out_path}
	cd ${out_path}
	mkdir -p ${prefix}
	cd ${prefix}
	WORKDIR=$(pwd)
	mkdir -p trimming assembly annotation NCBI logs

	## Trimming
	if [ ! -f ${WORKDIR}/trimming/trimming.done ];then
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\ttrimming the reads"
		fastp -i ${read1} -I ${read2} -o ${WORKDIR}/trimming/${base1} -O ${WORKDIR}/trimming/${base2} -w ${threads} &> ${WORKDIR}/logs/fastp.log
		touch ${WORKDIR}/trimming/trimming.done
	else 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\ttrimming was already done ... skipping"
	fi
	
	## de novo assembly
	if [ ! -f ${WORKDIR}/assembly/novoplasty.done ];then 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tperforming the de novo assembly"
		read1_novo=${WORKDIR}/trimming/${base1}
		read2_novo=${WORKDIR}/trimming/${base2}

		cat ${repo}/config_novo.txt |sed "s|WORKDIR|${WORKDIR}/assembly/|g"| sed "s|test|$prefix|" |sed "s|max_memory|${max_memory}|" | \
		sed "s|path_to_seed|${seed_path}|" |sed "s|range|$range|"|\
		sed "s|read1|$read1_novo|" | sed "s|read2|$read2_novo|" > ${prefix}_config.txt

		perl ${novo_path} -c ${WORKDIR}/${prefix}_config.txt &> ${WORKDIR}/logs/novoplasty.log

		touch ${WORKDIR}/assembly/novoplasty.done
	else 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tassembly was already done ... skipping"
	fi 

	## Standardizes the assembly
	if [ ! -f ${WORKDIR}/assembly/standardization.done ];then 
		${repo}/standardize_cpDNA.sh -d ${WORKDIR}/assembly/ -o ${WORKDIR}/assembly/${prefix}.plastome.fa -p ${prefix}
		touch ${WORKDIR}/assembly/standardization.done
	else 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tstandardization was already done ... skipping"
	fi 

	## Corrects non-ACTG characters in the assembly
	if [ ! -f ${WORKDIR}/assembly/check.done ];then 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tchecking for non-ATCG characters in the $prefix assembly"
		non_actg=$(cat ${WORKDIR}/assembly/${prefix}.plastome.fa|sed 1d|tr -d '\n' |tr -d 'ACTGactg')
		while [ ! -z "$non_actg" ];do 
			samtools faidx ${WORKDIR}/assembly/${prefix}.plastome.fa
			header=$(awk '{print $1}' ${WORKDIR}/assembly/${prefix}.plastome.fa.fai) 
			line1=$(cat ${WORKDIR}/assembly/${prefix}.plastome.fa|sed 1d|tr -d '\n' |sed 's/./&\n/g'|grep -n "[A-Za-z]"|grep -Ev "[ACTGactg]"|head -1|sed 's/:/\t/g')
			char=$(echo $line1|awk '{print $2}')
			if [ $(echo $line1|awk '{print $1}') -gt 100 ];then 
				loc=$(echo $line1|awk '{print "'$header':"$1-100"-"$1-1}') 
				seq=$(samtools faidx -n 200 ${WORKDIR}/assembly/${prefix}.plastome.fa "$loc"|sed 1d)
				new_seq=$(grep -o "$seq[A-Za-z]" ./sim/assembledGenome/Assembled_reads_sim_R*.fasta |cut -d ':' -f2|sort|uniq -c |sort -k1 -n -r |head -1|awk '{print $2}')
				sed -i "s/$seq$char/$new_seq/1" ${WORKDIR}/assembly/${prefix}.plastome.fa 
				non_actg=$(cat ${WORKDIR}/assembly/${prefix}.plastome.fa|sed 1d|tr -d '\n' |tr -d 'ACTGactg')
			else
				loc=$(echo $line1|awk '{print "'$header':"$1+1"-"$1+100}')
				seq=$(samtools faidx -n 200 ${WORKDIR}/assembly/${prefix}.plastome.fa "$loc"|sed 1d)
				new_seq=$(grep -o "[A-Za-z]$seq" ./sim/assembledGenome/Assembled_reads_sim_R*.fasta |cut -d ':' -f2|sort|uniq -c |sort -k1 -n -r |head -1|awk '{print $2}')
				sed -i "s/$char$seq/$new_seq/1" ${WORKDIR}/assembly/${prefix}.plastome.fa 
				non_actg=$(cat ${WORKDIR}/assembly/${prefix}.plastome.fa|sed 1d|tr -d '\n' |tr -d 'ACTGactg')
			fi 
		done 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\t$prefix Plastome assembly size: $(cat ${WORKDIR}/assembly/${prefix}.plastome.fa |sed 1d|tr -d '\n'|wc -c) bps"
		touch ${WORKDIR}/assembly/check.done
	else 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tchecking assembly quality was already done ... skipping"
	fi 

	## Annotation of the assembly
	if [ ! -f ${WORKDIR}/annotation/annotation.done ];then 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tannotating the assembly"
		cd ${WORKDIR}/annotation
		python3 ${repo}/AnnoPlast.py -f ${WORKDIR}/assembly/${prefix}.plastome.fa -g ${ref_gb_path} -o ./ -p ${prefix}
		touch ${WORKDIR}/annotation/annotation.done
	else 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tannotation was already done ... skipping"
	fi 

	## genbank to tbl format for NCBI submission
	if [ ! -f ${WORKDIR}/NCBI/tbl.done ];then 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tconverting genbank file to .tbl"
		cd ${WORKDIR}/NCBI/
		${repo}/gbf2tbl.pl ${WORKDIR}/annotation/${prefix}.gb
		mv ${WORKDIR}/annotation/${prefix}.tbl ${WORKDIR}/annotation/${prefix}.fsa . 
		touch ${WORKDIR}/NCBI/tbl.done
	else 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tgenbank to tbl conversion was already done ... skipping"
	fi 
	echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished plastaumatic on $prefix\n\n"
done 

