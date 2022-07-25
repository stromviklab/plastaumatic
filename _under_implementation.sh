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
repo=$(pwd)/$(dirname $0)

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







############### UNUSED
			## Selecting option1 or option2 based on the reference
			# for k in 1 2; do
			#     nucmer -p _tmp${k} -t ${threads} ${ref_fasta} ${WORKDIR}/assembledGenome/Option_${k}_${i}.standardized.fa
			#     if [ "$(show-tiling _tmp${k}.delta|sed 1d|awk '$5>95&&$7=="+" {print "yes"}')" == "yes" ];then
			#         cat ${WORKDIR}/assembledGenome/Option_${k}_${i}.standardized.fa|sed "s/Contig1/${prefix}/" > ${WORKDIR}/assembledGenome/${prefix}.plastome.fa
			#         echo "Selecting the Option_${k} based on the reference"
			#     fi
			#     rm -f _tmp${k}.delta
			# done

	# Annotation of the assembly
	# cd ${WORKDIR}/annotation
	# mkdir -p ref_plastome target_plstome out_plastome
	# cp ${ref_gb} ref_plastome
	# cp ${WORKDIR}/assembledGenome/${prefix}.plastome.fa target_plstome

	# #change the trna_length limitation in PGA.pl script (from 70 nt to 50 nt)
	# perl ${path_to_PGA} -r ref_plastome/ -t target_plstome/ -o out_plastome/

	# touch ${WORKDIR}/annotation/pga.done


	#check for ISC
	# ${repo}/check_internal_stops.sh ${WORKDIR}/annotation/out_plastome/${prefix}.plastome.mod.gb ${WORKDIR}/annotation/target_plstome/${prefix}.plastome.fa

	# Correction of ISC by correcting gene feature coordinates
	# dir=/scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly/PGA
	# ref=/scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly/genes.fa
	# cd ${dir}
	# mkdir corrections_ISC
	# cp ${dir}/out_plastome/*.gb corrections_ISC


	# correct gene features in ycf3 gene
	# for i in `cat ../ids`;do
	# echo $i
	# blastn -task blastn -subject ${dir}/target_plstome/${WORKDIR}/assembledGenome/${prefix}.plastome.mod.fa -query ${ref} -perc_identity 95 -evalue 10 -out ${i}.blast.ISC.out -outfmt '6 qseqid qstart qend sseqid sstart send length pident'

	# if [ $(cat ${i}.blast.ISC.out |grep "ycf3" |awk '$2=="1"&&$3==$7 {print $1}'|sort -u|wc -l) -eq 3 ];then
	# 	ycf3_3=$(cat ${i}.blast.ISC.out |grep "ycf3_3"|awk '$2=="1"&&$3==$7 {print $6".."$5}')
	# 	ycf3_2=$(cat ${i}.blast.ISC.out |grep "ycf3_2"|awk '$2=="1"&&$3==$7 {print $6".."$5}')
	# 	ycf3_1=$(cat ${i}.blast.ISC.out |grep "ycf3_1"|awk '$2=="1"&&$3==$7 {print $6".."$5}')
	# 	gene=$(echo "$ycf3_3..$ycf3_1"|awk -F '\\.\\.' '{print "     gene            complement("$1".."$4")"}')
	# 	cds=$(echo "     CDS             complement(join($ycf3_3,$ycf3_2,$ycf3_1))")
	# 	gene_old=$(grep -B1 "ycf3" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "gene ")
	# 	cds_old=$(grep -B1 "ycf3" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "CDS")
	# 	sed -i "s/$gene_old/$gene/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb
	# 	sed -i "s/$cds_old/$cds/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb
	# 	rm ${i}.blast.ISC.out
	# else
	# 	echo "Some of the CDS are missing alignments"
	# fi
	# done

	# correct gene feature rps12
	# for i in `cat ../ids`;do
	# 	echo $i
	# 	blastn -task blastn -subject ${dir}/target_plstome/${WORKDIR}/assembledGenome/${prefix}.plastome.mod.fa -query ${ref} \
	# 	-perc_identity 95 -evalue 10 -outfmt '6 qseqid qstart qend sseqid sstart send length \
	# 	pident sstrand' |grep "rps12" |awk '$2=="1"&&$3==$7 {print}' > ${i}.blast.ISC.out

	# 	if [ $(cat ${i}.blast.ISC.out |cut -f1|sort -u|wc -l) -eq 5 ];then

	#     #rps12 first gene
	#     rps12_3=$(cat ${i}.blast.ISC.out |grep "rps12_1_3"|grep "minus"|awk '{print $6".."$5}')
	#     rps12_2=$(cat ${i}.blast.ISC.out |grep "rps12_1_2"|grep "minus"|awk '{print $6".."$5}')
	#     rps12_1=$(cat ${i}.blast.ISC.out |grep "rps12_12_1"|grep "minus"|awk '{print $6".."$5}')
	#     gene=$(echo "$rps12_3..$rps12_2"|awk -F '\\.\\.' '{print "     gene            complement(join("$1".."$4",""'$rps12_1'""))"}')
	#     cds=$(echo "     CDS             complement(join($rps12_3,$rps12_2,$rps12_1))")
	#     gene_old=$(grep -B1 "rps12" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "gene "|head -1)
	#     cds_old=$(grep -B1 "rps12" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "CDS"|head -1)
	#     sed -i "s/$gene_old/$gene/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb
	#     sed -i "s/$cds_old/$cds/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb

	#     #rps12 second gene
	#     rps12_3=$(cat ${i}.blast.ISC.out |grep "rps12_2_3"|grep "plus"|awk '{print $5".."$6}')
	#     rps12_2=$(cat ${i}.blast.ISC.out |grep "rps12_2_2"|grep "plus"|awk '{print $5".."$6}')
	#     rps12_1=$(cat ${i}.blast.ISC.out |grep "rps12_12_1"|grep "minus"|awk '{print $6".."$5}')
	#     gene=$(echo "$rps12_2..$rps12_3"|awk -F '\\.\\.' '{print "     gene            join(complement(""'$rps12_1'""),"$1".."$4")"}')
	#     cds=$(echo "     CDS             join(complement($rps12_1),$rps12_2,$rps12_3)")
	#     gene_old=$(grep -B1 "rps12" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "gene "|head -2|tail -1)
	#     cds_old=$(grep -B1 "rps12" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "CDS"|head -2|tail -1)
	#     sed -i "s/$gene_old/$gene/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb
	#     sed -i "s/$cds_old/$cds/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb

	#     del=$(grep -B1 "rps12" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "gene "|tail -1)
	#     sed -i "/$del/,+6d" ${dir}/corrections_ISC/${i}.plastome.mod.gb
	#     rm ${i}.blast.ISC.out
	# 	else
	# 		echo "Some of the CDS are missing alignments"
	# 	fi
	# done

	# # add exceptions and notes
	# for i in `cat ../ids`;do
	# 	echo $i
	# 	sed -i 's/gene=\"rps12\"/gene=\"rps12\"\n                     \/trans_splicing/g' ${dir}/corrections_ISC/${i}.plastome.mod.gb
	# 	sed -i 's/gene=\"rps19\"/gene=\"rps19\"\n                     \/note=\"GTG start codon\"/g' ${dir}/corrections_ISC/${i}.plastome.mod.gb
	# 	sed -i 's/gene=\"ndhD\"/gene=\"ndhD\"\n                     \/exception=\"RNA editing\"/g' ${dir}/corrections_ISC/${i}.plastome.mod.gb
	# 	sed -i 's/gene=\"psbL\"/gene=\"psbL\"\n                     \/exception=\"RNA editing\"/g' ${dir}/corrections_ISC/${i}.plastome.mod.gb
	# done


	# # Modify the headers in fasta and genbank for NCBI submission
	# # get the attributes file with organims name, bioproject and biosample
	# cd /scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly
	# mkdir NCBI_sub && cd NCBI_sub
	# mkdir fasta_plastome tbl_plastome
	# grep ">" ../PGA/target_plstome/*.fa|cut -d ':' -f2 |sort > a
	# cat attributes.txt |awk -F '\t' '{print ">"$1" [organism="$2"] [Bioproject="$3"] [Biosample="$4"]"}' OFS='\t' |sort > attributes.mod.txt
	# paste a attributes.mod.txt > b
	# mv b attributes.mod.txt
	# rm a
	# cat  attributes.mod.txt |while read i;do
	# 	id=$(echo "$i"|awk -F '\t' '{print $1}'|tr -d '>')
	# 	old=$(echo "$i"|awk -F '\t' '{print $1}')
	# 	new=$(echo "$i"|awk -F '\t' '{print $2}')
	# 	cat /scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly/PGA/target_plstome/${id}*plastome.mod.fa|\
	# 	sed "s/$old/$new/1" > fasta_plastome/${id}_plastome_assembly.fa
	# 	old2=$(echo "$id.plastome.mod")
	# 	new2=$(echo "$i"|awk -F '\t' '{print $2}'|cut -d ' ' -f1|tr -d '>')
	# 	cat /scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly/PGA/tbl/${id}*plastome.mod.tbl|\
	# 	sed "s/$old2/$new2/1" > tbl_plastome/${id}_plastome_assembly.tbl
	# done

	# cat fasta_plastome/*.fa  > ICP_plastomes_NCBI_submission.fa
	# cat tbl_plastome/*.tbl  > ICP_plastomes_NCBI_submission.tbl
