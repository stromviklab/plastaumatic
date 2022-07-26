#!/bin/bash

## Exits if a command fails
set -eo pipefail

## Checks if the script is run properly
if [ $# -ne 1 ];then
	echo "Usage:
	job_input.sh test.cfg"
	exit 1
fi

## Checks if an input config file exists
if [ ! -f $1 ];then
	echo -e "ERROR: \t Config file does not exist"
	exit 1
fi

source $1

## Checks if the input files are present
if [ ! -f $seed ];then
	echo -e "ERROR:\t Seed file does not exist, check again" && exit 1
elif [ ! -f $ref_gb ];then 
	echo -e "ERROR:\t Reference GenBank file does not exist, check again" && exit 1
elif [ ! -d $out ];then 
	echo -e "ERROR:\t Output directory does not exist, check again" && exit 1
elif [ ! -f $path_to_trimmomatic ];then 
	echo -e "ERROR:\t Trimmomatic jar file does not exist, check again" && exit 1
elif [ ! -f $adapters ];then 
	echo -e "ERROR:\t Adapter fasta file does not exist, check again" && exit 1
elif [ ! -f $path_to_novoplasty ];then 
	echo -e "ERROR:\t NOVOPlasty executable does not exist, check again" && exit 1
elif [ ! -f $fof ];then 
	echo -e "ERROR:\t fof file does not exist, check again" && exit 1
fi

## Sets memory and cpus available
# threads=$(nproc --all)
# max_memory=$(awk '/MemAvailable/ {printf "%.0f\n", $2/1024/1024*0.9}' /proc/meminfo)
threads=1
max_memory=30
repo=$(dirname $0)

for line in `cat $fof`;do 
	prefix=$(echo $line|cut -d ',' -f1)
	read1=$(echo $line|cut -d ',' -f2)
	read2=$(echo $line|cut -d ',' -f3)

  if [ ! -f $read1 ] || [ ! -f $read2 ];then
		echo -e "ERROR: \t Read files does not exist, check again" && exit 1
	fi  

	## Setting the directories
	echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tsetting up the directories" 
	cd ${out}
	mkdir -p ${prefix}
	cd ${prefix}
	WORKDIR=$(pwd)
	mkdir -p reads trimmedReads assembledGenome annotation NCBI logs

	## decompress gzipped reads
	if [ ! -f ${WORKDIR}/reads/reads.done ];then 
		if [ $(echo "${read1##*.}") == "gz" ];then
			echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tdecompressing the read files"
			pigz -d -p ${threads} -k -c ${read1} > ${WORKDIR}/reads/${prefix}_1.fq
			pigz -d -p ${threads} -k -c ${read2} > ${WORKDIR}/reads/${prefix}_2.fq
		else
			ln -sf ${read1} ${WORKDIR}/reads/${prefix}_1.fq
			ln -sf ${read2} ${WORKDIR}/reads/${prefix}_2.fq
		fi
		touch ${WORKDIR}/reads/reads.done
	else
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tdecompressing was already done ... skipping"
	fi

	## Trimming
	if [ ! -f ${WORKDIR}/trimmedReads/trimming.done ];then
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\ttrimming the reads"
		java -jar ${path_to_trimmomatic} PE -threads ${threads} ${WORKDIR}/reads/${prefix}_1.fq \
		${WORKDIR}/reads/${prefix}_2.fq ${WORKDIR}/trimmedReads/${prefix}_1P.fq ${WORKDIR}/trimmedReads/${prefix}_1U.fq \
		${WORKDIR}/trimmedReads/${prefix}_2P.fq ${WORKDIR}/trimmedReads/${prefix}_2U.fq \
		ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60 &> ${WORKDIR}/logs/trimmomatic.log
		touch ${WORKDIR}/trimmedReads/trimming.done
	else 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\ttrimming was already done ... skipping"
	fi
	

	## de novo assembly
	if [ ! -f ${WORKDIR}/assembledGenome/novoplasty.done ];then 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tperforming the de novo assembly"
		read1_novo=${WORKDIR}/trimmedReads/${prefix}_1P.fq
		read2_novo=${WORKDIR}/trimmedReads/${prefix}_2P.fq

		cat ${repo}/config_novo.txt |sed "s|WORKDIR|${WORKDIR}/assembledGenome/|g"| sed "s|test|$prefix|" |sed "s|max_memory|${max_memory}|" | \
		sed "s|path_to_seed|${seed}|" |sed "s|range|$range|"|\
		sed "s|read1|$read1_novo|" | sed "s|read2|$read2_novo|" > ${prefix}_config.txt

		perl ${path_to_novoplasty} -c ${WORKDIR}/${prefix}_config.txt &> ${WORKDIR}/logs/novoplasty.log

		touch ${WORKDIR}/assembledGenome/novoplasty.done
	else 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tassembly was already done ... skipping"
	fi 

	## Standardizes the assembly
	if [ ! -f ${WORKDIR}/assembledGenome/standardization.done ];then 
		${repo}/standardize_cpDNA.sh -d ${WORKDIR}/assembledGenome/ -o ${WORKDIR}/assembledGenome/${prefix}.plastome.fa -p ${prefix}
		touch ${WORKDIR}/assembledGenome/standardization.done
	else 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tstandardization was already done ... skipping"
	fi 

	## Corrects non-ACTG characters in the assembly
	if [ ! -f ${WORKDIR}/assembledGenome/check.done ];then 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tchecking for non-ATCG characters in the $prefix assembly"
		non_actg=$(cat ${WORKDIR}/assembledGenome/${prefix}.plastome.fa|sed 1d|tr -d '\n' |tr -d 'ACTGactg')
		while [ ! -z "$non_actg" ];do 
			samtools faidx ${WORKDIR}/assembledGenome/${prefix}.plastome.fa
			header=$(awk '{print $1}' ${WORKDIR}/assembledGenome/${prefix}.plastome.fa.fai) 
			line1=$(cat ${WORKDIR}/assembledGenome/${prefix}.plastome.fa|sed 1d|tr -d '\n' |sed 's/./&\n/g'|grep -n "[A-Za-z]"|grep -Ev "[ACTGactg]"|head -1|sed 's/:/\t/g')
			char=$(echo $line1|awk '{print $2}')
			if [ $(echo $line1|awk '{print $1}') -gt 100 ];then 
				loc=$(echo $line1|awk '{print "'$header':"$1-100"-"$1-1}') 
				seq=$(samtools faidx -n 200 ${WORKDIR}/assembledGenome/${prefix}.plastome.fa "$loc"|sed 1d)
				new_seq=$(grep -o "$seq[A-Za-z]" ./sim/assembledGenome/Assembled_reads_sim_R*.fasta |cut -d ':' -f2|sort|uniq -c |sort -k1 -n -r |head -1|awk '{print $2}')
				sed -i "s/$seq$char/$new_seq/1" ${WORKDIR}/assembledGenome/${prefix}.plastome.fa 
				non_actg=$(cat ${WORKDIR}/assembledGenome/${prefix}.plastome.fa|sed 1d|tr -d '\n' |tr -d 'ACTGactg')
			else
				loc=$(echo $line1|awk '{print "'$header':"$1+1"-"$1+100}')
				seq=$(samtools faidx -n 200 ${WORKDIR}/assembledGenome/${prefix}.plastome.fa "$loc"|sed 1d)
				new_seq=$(grep -o "[A-Za-z]$seq" ./sim/assembledGenome/Assembled_reads_sim_R*.fasta |cut -d ':' -f2|sort|uniq -c |sort -k1 -n -r |head -1|awk '{print $2}')
				sed -i "s/$char$seq/$new_seq/1" ${WORKDIR}/assembledGenome/${prefix}.plastome.fa 
				non_actg=$(cat ${WORKDIR}/assembledGenome/${prefix}.plastome.fa|sed 1d|tr -d '\n' |tr -d 'ACTGactg')
			fi 
		done 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\t$prefix Plastome assembly size: $(cat ${WORKDIR}/assembledGenome/${prefix}.plastome.fa |sed 1d|tr -d '\n'|wc -c) bps"
		touch ${WORKDIR}/assembledGenome/check.done
	else 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tchecking assembly quality was already done ... skipping"
	fi 

	## Annotation of the assembly
	if [ ! -f ${WORKDIR}/annotation/annotation.done ];then 
		echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tannotating the assembly"
		cd ${WORKDIR}/annotation
		python3 ${repo}/AnnoPlast.py -f ${WORKDIR}/assembledGenome/${prefix}.plastome.fa -g ${ref_gb} -o ./ -p ${prefix}
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
done 

