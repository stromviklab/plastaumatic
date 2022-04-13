#!/bin/bash

# Checks if the script is run properly
if [ $# -ne 1 ];then
    echo "Usage: 
        job_input.sh test.cfg"
    exit 1
fi

# Checks if an input config file exists
if [ ! -f $1 ];then 
  echo "config file does not exist"
else 
    source $1
    cd ${out}
    mkdir -p ${prefix} 
    cd ${prefix}
    WORKDIR=$(pwd)
    mkdir -p reads trimmedReads assembly annotation NCBI logs 
    if [ $(echo "${read1##*.}") == "gz" ];then 
        pigz -d -p ${threads_available} -k -c ${read1} > ${WORKDIR}/reads/${prefix}_1.fq 
        pigz -d -p ${threads_available} -k -c ${read2} > ${WORKDIR}/reads/${prefix}_2.fq
    else 
        ln -sf ${read1} ${WORKDIR}/reads/${prefix}_1.fq
        ln -sf ${read2} ${WORKDIR}/reads/${prefix}_2.fq
    fi
    adaptor="$(dirname "$path_to_trimmomatic")"/adapters/
    if [ -d $adaptor ];then 
        java -jar ${path_to_trimmomatic} PE -threads ${threads_available} ${WORKDIR}/reads/${prefix}_1.fq \
        ${WORKDIR}/reads/${prefix}_2.fq ${WORKDIR}/trimmedReads/${prefix}_1P.fq ${WORKDIR}/trimmedReads/${prefix}_1U.fq \
        ${WORKDIR}/trimmedReads/${prefix}_2P.fq ${WORKDIR}/trimmedReads/${prefix}_2U.fq \
        ILLUMINACLIP:${adaptor}/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60 &> ${WORKDIR}/logs/trimmomatic.log 
    else 
        echo "ERROR: can't find the trimmomatic adapters folder"
    fi 

    read1_novo=${WORKDIR}/trimmedReads/${prefix}_1P.fq
    read2_novo=${WORKDIR}/trimmedReads/${prefix}_2P.fq

    cat ${path_to_repo}/config_novo.txt |sed "s|WORKDIR|${WORKDIR}/assembly|g"| sed "s|test|$prefix|" |sed "s|max_memory|${max_memory}|" | \
    sed "s|path_to_seed|${seed}|" |sed "s|path_to_reference|${ref_fasta}|" |sed "s|range|$range|"|\
    sed "s|read1|$read1_novo|" | sed "s|read2|$read2_novo|" > ${prefix}_config.txt


    perl ${path_to_novoplasty} -c ${WORKDIR}/${prefix}_config.txt &> ${WORKDIR}/logs/novoplasty.log 

    if [ -f ${WORKDIR}/assembly/Option_1_${prefix}.fasta ];then 
        ${path_to_repo}/standardize_cpDNA.sh -i ${WORKDIR}/assembly/Option_1_${prefix}.fasta
        ${path_to_repo}/standardize_cpDNA.sh -i ${WORKDIR}/assembly/Option_2_${prefix}.fasta
        for k in 1 2; do
            nucmer -p _tmp${k} -t ${threads_available} ${ref_fasta} ${WORKDIR}/assembly/Option_${k}_${i}.standardized.fa
            if [ "$(show-tiling _tmp${k}.delta|sed 1d|awk '$5>95&&$7=="+" {print "yes"}')" == "yes" ];then 
                cat ${WORKDIR}/assembly/Option_${k}_${i}.standardized.fa|sed "s/Contig1/${prefix}/" > ${WORKDIR}/assembly/${prefix}.plastome.fa 
                echo "Selecting the Option_${k} based on the reference"
            fi 
            rm -f _tmp${k}.delta 
        done 
    else 
        if [ -f ${WORKDIR}/assembly/Circularized_assembly_1_${prefix}.fasta ];then
            ${path_to_repo}/standardize_cpDNA.sh -i ${WORKDIR}/assembly/Circularized_assembly_1_${prefix}.fasta
            cat ${WORKDIR}/assembly/Circularized_assembly_1_${i}.standardized.fa|sed "s/Contig1/${prefix}/" > ${WORKDIR}/assembly/${prefix}.plastome.fa 
            echo "Since only one assembly is avaialble, selecting that only"
        fi
    fi

    cp ${WORKDIR}/assembly/${prefix}.plastome.fa ${WORKDIR}/assembly/${prefix}.plastome.mod.fa
    sed -i ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' ${WORKDIR}/assembly/${prefix}.plastome.mod.fa
    samtools faidx ${WORKDIR}/assembly/${prefix}.plastome.mod.fa 
    lines=$(cat ${WORKDIR}/assembly/${prefix}.plastome.mod.fa |sed 1d|tr -d 'ACTG'|tr -d 'actg'|grep -v "^$"|wc -l)
    if [ $(echo $lines) -gt 0 ];then 
        cat ${WORKDIR}/assembly/${prefix}.plastome.mod.fa |sed 1d|tr -d 'ACTG'|tr -d 'actg'|grep -v "^$"|grep -o . |sort -u|\
        while read char;do 
            seqkit locate  --pattern "$char" ${WORKDIR}/assembly/${prefix}.plastome.mod.fa|sed 1d|\
            while read loc;do 
                loc1=$(echo $loc|awk '{print $1":"$5-50"-"$5-1}')
                loc2=$(echo $loc|awk '{print $1":"$5+1"-"$6+50}')
                seq1=$(samtools faidx -n 200 ${WORKDIR}/assembly/${prefix}.plastome.mod.fa "$loc1"|sed 1d)
                seq2=$(samtools faidx -n 200 ${WORKDIR}/assembly/${prefix}.plastome.mod.fa "$loc2"|sed 1d)
                new_seq=$(grep -o "$seq1[A-Za-z]$seq2" ${WORKDIR}/trimmedReads/${prefix}_1P.fq |sort|uniq -c |sort -k1 -n -r |head -1|awk '{print $2}')
                sed -i "s/$seq1$char$seq2/$new_seq/g" ${WORKDIR}/assembly/${prefix}.plastome.mod.fa
            done 
        done 
        lines2=$(cat ${WORKDIR}/assembly/${prefix}.plastome.mod.fa |sed 1d|tr -d 'ACTG'|tr -d 'actg'|grep -v "^$"|wc -l)
        size1=$(cat ${WORKDIR}/assembly/${prefix}.plastome.fa |sed 1d|tr -d '\n'|wc -c)
        size2=$(cat ${WORKDIR}/assembly/${prefix}.plastome.mod.fa |sed 1d|tr -d '\n'|wc -c)
        if [[ $(echo $lines2) -eq 0 && $size1 -eq $size2 ]];then 
            echo -e "$prefix:\n  non-ATCG characters were corrected\n  Assembly size: $size2"
        else 
            echo "Something went wrong"
        fi 
    else 
        echo -e "$prefix:\n  No non-ATCG characters are found\n  Assembly size: $(cat ${WORKDIR}/assembly/${prefix}.plastome.fa |sed 1d|tr -d '\n'|wc -c)"
    fi 
    #annotate
    cd ${WORKDIR}/annotation
    mkdir -p ref_plastome target_plstome out_plastome
    cp ${ref_gb} ref_plastome 
    cp ${WORKDIR}/assembly/${prefix}.plastome.mod.fa target_plstome

    #change the trna_length limitation in PGA.pl script (from 70 nt to 50 nt)
    perl ${path_to_repo}/PGA.pl -r ref_plastome/ -t target_plstome/ -o out_plastome/

    #check for ISC
    ${path_to_repo}/check_internal_stops.sh ${WORKDIR}/annotation/out_plastome/${prefix}.plastome.mod.gb ${WORKDIR}/annotation/target_plstome/${prefix}.plastome.mod.fa

# Correction of ISC by correcting gene feature coordinates 
dir=/scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly/PGA
ref=/scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly/genes.fa
cd ${dir}
mkdir corrections_ISC
cp ${dir}/out_plastome/*.gb corrections_ISC
# correct gene features in ycf3 gene
for i in `cat ../ids`;do 
    echo $i
    blastn -task blastn -subject ${dir}/target_plstome/${WORKDIR}/assembly/${prefix}.plastome.mod.fa -query ${ref} -perc_identity 95 -evalue 10 -out ${i}.blast.ISC.out -outfmt '6 qseqid qstart qend sseqid sstart send length pident'

    if [ $(cat ${i}.blast.ISC.out |grep "ycf3" |awk '$2=="1"&&$3==$7 {print $1}'|sort -u|wc -l) -eq 3 ];then 
        ycf3_3=$(cat ${i}.blast.ISC.out |grep "ycf3_3"|awk '$2=="1"&&$3==$7 {print $6".."$5}')
        ycf3_2=$(cat ${i}.blast.ISC.out |grep "ycf3_2"|awk '$2=="1"&&$3==$7 {print $6".."$5}')
        ycf3_1=$(cat ${i}.blast.ISC.out |grep "ycf3_1"|awk '$2=="1"&&$3==$7 {print $6".."$5}')
        gene=$(echo "$ycf3_3..$ycf3_1"|awk -F '\\.\\.' '{print "     gene            complement("$1".."$4")"}')
        cds=$(echo "     CDS             complement(join($ycf3_3,$ycf3_2,$ycf3_1))")
        gene_old=$(grep -B1 "ycf3" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "gene ")
        cds_old=$(grep -B1 "ycf3" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "CDS")
        sed -i "s/$gene_old/$gene/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb
        sed -i "s/$cds_old/$cds/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb
        rm ${i}.blast.ISC.out
    else
        echo "Some of the CDS are missing alignments"
    fi 
done 

# correct gene feature rps12
for i in `cat ../ids`;do 
    echo $i
    blastn -task blastn -subject ${dir}/target_plstome/${WORKDIR}/assembly/${prefix}.plastome.mod.fa -query ${ref} \
    -perc_identity 95 -evalue 10 -outfmt '6 qseqid qstart qend sseqid sstart send length \
    pident sstrand' |grep "rps12" |awk '$2=="1"&&$3==$7 {print}' > ${i}.blast.ISC.out

    if [ $(cat ${i}.blast.ISC.out |cut -f1|sort -u|wc -l) -eq 5 ];then 

        #rps12 first gene
        rps12_3=$(cat ${i}.blast.ISC.out |grep "rps12_1_3"|grep "minus"|awk '{print $6".."$5}')
        rps12_2=$(cat ${i}.blast.ISC.out |grep "rps12_1_2"|grep "minus"|awk '{print $6".."$5}')
        rps12_1=$(cat ${i}.blast.ISC.out |grep "rps12_12_1"|grep "minus"|awk '{print $6".."$5}')
        gene=$(echo "$rps12_3..$rps12_2"|awk -F '\\.\\.' '{print "     gene            complement(join("$1".."$4",""'$rps12_1'""))"}')
        cds=$(echo "     CDS             complement(join($rps12_3,$rps12_2,$rps12_1))")
        gene_old=$(grep -B1 "rps12" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "gene "|head -1)
        cds_old=$(grep -B1 "rps12" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "CDS"|head -1)
        sed -i "s/$gene_old/$gene/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb
        sed -i "s/$cds_old/$cds/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb

        #rps12 second gene 
        rps12_3=$(cat ${i}.blast.ISC.out |grep "rps12_2_3"|grep "plus"|awk '{print $5".."$6}')
        rps12_2=$(cat ${i}.blast.ISC.out |grep "rps12_2_2"|grep "plus"|awk '{print $5".."$6}')
        rps12_1=$(cat ${i}.blast.ISC.out |grep "rps12_12_1"|grep "minus"|awk '{print $6".."$5}')
        gene=$(echo "$rps12_2..$rps12_3"|awk -F '\\.\\.' '{print "     gene            join(complement(""'$rps12_1'""),"$1".."$4")"}')
        cds=$(echo "     CDS             join(complement($rps12_1),$rps12_2,$rps12_3)")
        gene_old=$(grep -B1 "rps12" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "gene "|head -2|tail -1)
        cds_old=$(grep -B1 "rps12" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "CDS"|head -2|tail -1)
        sed -i "s/$gene_old/$gene/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb
        sed -i "s/$cds_old/$cds/1" ${dir}/corrections_ISC/${i}.plastome.mod.gb

        del=$(grep -B1 "rps12" ${dir}/corrections_ISC/${i}.plastome.mod.gb |grep -wF "gene "|tail -1)
        sed -i "/$del/,+6d" ${dir}/corrections_ISC/${i}.plastome.mod.gb
        rm ${i}.blast.ISC.out
    else
        echo "Some of the CDS are missing alignments"
    fi 
done 

# add exceptions and notes
for i in `cat ../ids`;do 
    echo $i
    sed -i 's/gene=\"rps12\"/gene=\"rps12\"\n                     \/trans_splicing/g' ${dir}/corrections_ISC/${i}.plastome.mod.gb
    sed -i 's/gene=\"rps19\"/gene=\"rps19\"\n                     \/note=\"GTG start codon\"/g' ${dir}/corrections_ISC/${i}.plastome.mod.gb
    sed -i 's/gene=\"ndhD\"/gene=\"ndhD\"\n                     \/exception=\"RNA editing\"/g' ${dir}/corrections_ISC/${i}.plastome.mod.gb
    sed -i 's/gene=\"psbL\"/gene=\"psbL\"\n                     \/exception=\"RNA editing\"/g' ${dir}/corrections_ISC/${i}.plastome.mod.gb
done 

# Use gbf2tbl.pl script to convert to tbl format and submit to NCBI
# wget ftp://ftp.ncbi.nlm.nih.gov//toolbox/ncbi_tools/converters/scripts/gbf2tbl.pl
cd /scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly/PGA/corrections_ISC
scr=/scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly/PGA/gbf2tbl.pl
for i in `cat /scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly/ids`;do 
    echo $i
    ${scr} $i.plastome.mod.gb
done 
rm *.fsa 
cd ../
mkdir tbl 
mv corrections_ISC/*.tbl tbl/ 

# Modify the headers in fasta and genbank for NCBI submission 
# get the attributes file with organims name, bioproject and biosample
cd /scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly
mkdir NCBI_sub && cd NCBI_sub
mkdir fasta_plastome tbl_plastome
grep ">" ../PGA/target_plstome/*.fa|cut -d ':' -f2 |sort > a
cat attributes.txt |awk -F '\t' '{print ">"$1" [organism="$2"] [Bioproject="$3"] [Biosample="$4"]"}' OFS='\t' |sort > attributes.mod.txt 
paste a attributes.mod.txt > b
mv b attributes.mod.txt 
rm a 
cat  attributes.mod.txt |while read i;do 
    id=$(echo "$i"|awk -F '\t' '{print $1}'|tr -d '>')
    old=$(echo "$i"|awk -F '\t' '{print $1}')
    new=$(echo "$i"|awk -F '\t' '{print $2}')
    cat /scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly/PGA/target_plstome/${id}*plastome.mod.fa|\
    sed "s/$old/$new/1" > fasta_plastome/${id}_plastome_assembly.fa
    old2=$(echo "$id.plastome.mod")
    new2=$(echo "$i"|awk -F '\t' '{print $2}'|cut -d ' ' -f1|tr -d '>')
    cat /scratch/m/mstrom/reddy7/ICP-2021/organelle_assembly/PGA/tbl/${id}*plastome.mod.tbl|\
    sed "s/$old2/$new2/1" > tbl_plastome/${id}_plastome_assembly.tbl 
 done 

cat fasta_plastome/*.fa  > ICP_plastomes_NCBI_submission.fa
cat tbl_plastome/*.tbl  > ICP_plastomes_NCBI_submission.tbl 

fi 
