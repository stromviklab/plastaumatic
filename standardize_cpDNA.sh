#!/bin/bash

while getopts 'd:o:p:' options
do
  case $options in
    d) dir=$OPTARG ;;
    o) out=$OPTARG ;;
    p) prefix=$OPTARG ;;
  esac
done

if [ -f ${dir}/Circularized_assembly_1_${prefix}.fasta ];then
  fasta=${dir}/Circularized_assembly_1_${prefix}.fasta
elif [ -f ${dir}/Option_1_${prefix}.fasta ];then 
  fasta=${dir}/Option_1_${prefix}.fasta
elif [ -f ${dir}/Option_2_${prefix}.fasta ];then 
  fasta=${dir}/Option_2_${prefix}.fasta
else 
  echo "ERROR: NOVOPlasty could not return a complete circular assembly, assemble manually"
  exit 1
fi 

echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tstandardizing the input fasta file"

## index the fasta input and get the sequence length
samtools faidx ${fasta}
seq_len=$(cut -f2 ${fasta}.fai)

mkdir -p _tmp_$(basename $fasta) && cd _tmp_$(basename $fasta)
## get the header 
head -1 ${fasta} > header

## blast the fasta to itself to get the repeat alignments
## ideally, the IR sequences in plastome should be identical, and hence a perc_identity of 99 is used
## self match is removed and any match less than 1000 bp length is removed
## based on test runs, the length can be reduced down to 100 bp (not sure how it preforms with 100 though)
blastn -query ${fasta} -subject ${fasta} -perc_identity 99 -evalue 0.00001 -outfmt '6 qseqid qstart qend sseqid sstart send length pident'|awk '$2!=$5&&$3!=$6 {print}'|awk '$7>1000 {print}' > _tmp.blast.out

## if the fasta is split at inverted repeats, it would give more than 2 alignments (4 to be precise)
if [ $(grep -c "^" _tmp.blast.out) -gt 2 ];then 
    ir_tmp=$(cat _tmp.blast.out|awk '$3=='$seq_len' {print $6-1}')
    sc1=$(cat _tmp.blast.out|awk '$5=='$ir_tmp' {print $6-$3}')
    sc2=$(cat _tmp.blast.out|awk '$3=='$seq_len' {print $2-$5}')
    if [ $sc2 -gt $sc1 ];then 
        ## do this if the fasta is split at IRa region 
        lsc=$(cat _tmp.blast.out|awk '$3=='$seq_len' {print $1":"$5+1"-"$2-1}')
        ira_start=$(cat _tmp.blast.out|awk '$3=='$seq_len' {print $1":"$2"-"$3}')
        ira_end=$(cat _tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$2"-"$3}') 
        ssc=$(cat _tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$3+1"-"$6-1}')
        irb_start=$(cat _tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$6"-"$5}') 
        irb_end=$(cat _tmp.blast.out|awk '$3=='$seq_len' {print $1":"$6"-"$5}')
    else 
        ## do this if the fasta is split at IRb region 
        lsc=$(cat _tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$3+1"-"$6-1}')
        ira_start=$(cat _tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$6"-"$5}') 
        ira_end=$(cat _tmp.blast.out|awk '$3=='$seq_len' {print $1":"$6"-"$5}')
        ssc=$(cat _tmp.blast.out|awk '$3=='$seq_len' {print $1":"$5+1"-"$2-1}')
        irb_start=$(cat _tmp.blast.out|awk '$3=='$seq_len' {print $1":"$2"-"$3}')
        irb_end=$(cat _tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$2"-"$3}') 
    fi 
    samtools faidx ${fasta} "$irb_start" |sed "s/$irb_start/IRb/" > irbs 
    samtools faidx ${fasta} "$irb_end" |sed 1d > irbe  
    samtools faidx ${fasta} "$ira_start" |sed "s/$ira_start/IRa/" > iras 
    samtools faidx ${fasta} "$ira_end" |sed 1d > irae  
    cat irbs irbe > IRb.fa 
    cat iras irae > IRa.fa 
    ## double checking the IR alignments 
    if [ "$(blastn -query IRa.fa -subject IRb.fa -perc_identity 99 -evalue 0.00001 -outfmt '6 qseqid qstart qend sseqid sstart send length pident' |awk 'NR==1&&$2==$6&&$3==$5 {print "Success"}')" == "Success" ] ;then 
      samtools faidx ${fasta} "$lsc" |sed 1d > LSC.fa 
      samtools faidx ${fasta} "$ssc" |sed 1d > SSC.fa 
      cat IRa.fa | sed 1d > _tmp.IRa.fa 
      cat IRb.fa | sed 1d > _tmp.IRb.fa
      cat header LSC.fa _tmp.IRa.fa SSC.fa _tmp.IRb.fa | sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' | sed "s/^>.*$/>${prefix}/" > ${out}
      rm header LSC.fa SSC.fa _tmp.IRa.fa _tmp.IRb.fa IRa.fa IRb.fa _tmp.blast.out irbs irbe iras irae ${fasta}.fai 
      echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
    fi

else 
    ## if the end or start the sequence is there in the alignemnts, that means the IR locations needs to be adjusted
    if [ "$(cat _tmp.blast.out| awk '$3=='$seq_len'||$2==1 {print "yes"}')" == "yes" ];then 
      if [ "$(cat _tmp.blast.out| awk '$3=='$seq_len' {print "yes"}')" == "yes" ];then
        if [ $(awk '$3=='$seq_len' {print $6}' _tmp.blast.out) -gt 50000 ];then
          sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' ${fasta} | sed "s/^>.*$/>${prefix}/" > ${out}
          rm _tmp.blast.out ${fasta}.fai header
          echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
        else 
          ira=$(cat _tmp.blast.out |awk '$3=='$seq_len' {print $1":"$2"-"$3}')
          irb=$(cat _tmp.blast.out |awk '$3=='$seq_len' {print $1":"$6"-"$5}')
          lsc=$(cat _tmp.blast.out |awk '$3=='$seq_len' {print $1":"$5+1"-"$2-1}')  
          ssc=$(cat _tmp.blast.out |awk '$3=='$seq_len' {print $1":1-"$6-1}')  

          samtools faidx ${fasta} "$irb" |sed 1d > IRb.fa
          samtools faidx ${fasta} "$ira" |sed 1d > IRa.fa 
          samtools faidx ${fasta} "$ssc" |sed 1d > ssc.fa 
          samtools faidx ${fasta} "$lsc" |sed 1d > lsc.fa 

          cat header lsc.fa IRa.fa ssc.fa IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' | sed "s/^>.*$/>${prefix}/" > ${out}
          rm _tmp.blast.out lsc.fa ssc.fa IRa.fa IRb.fa header ${fasta}.fai
          echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
        fi 
      fi 
      if [ "$(cat _tmp.blast.out| awk '$2==1 {print "yes"}')" == "yes" ];then
        if [ $(cat _tmp.blast.out |awk '$2==1 {print '$seq_len'-$5}') -gt 50000 ];then 
          ira=$(cat _tmp.blast.out |awk '$2==1 {print $1":"$2"-"$3}')
          irb=$(cat _tmp.blast.out |awk '$2==1 {print $1":"$6"-"$5}')
          lsc=$(cat _tmp.blast.out |awk '$2==1 {print $1":"$5+1"-"'$seq_len'}')  
          ssc=$(cat _tmp.blast.out |awk '$2==1 {print $1":"$3+1"-"$6-1}')  

          samtools faidx ${fasta} "$irb" |sed 1d > IRb.fa
          samtools faidx ${fasta} "$ira" |sed 1d > IRa.fa 
          samtools faidx ${fasta} "$ssc" |sed 1d > ssc.fa 
          samtools faidx ${fasta} "$lsc" |sed 1d > lsc.fa 

          cat header lsc.fa IRa.fa ssc.fa IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'| sed "s/^>.*$/>${prefix}/" > ${out}
          rm _tmp.blast.out lsc.fa ssc.fa IRa.fa IRb.fa header ${fasta}.fai
          echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
        else 
          irb=$(cat _tmp.blast.out |awk '$2==1 {print $1":"$2"-"$3}')
          ira=$(cat _tmp.blast.out |awk '$2==1 {print $1":"$6"-"$5}')
          ssc=$(cat _tmp.blast.out |awk '$2==1 {print $1":"$5+1"-"'$seq_len'}')  
          lsc=$(cat _tmp.blast.out |awk '$2==1 {print $1":"$3+1"-"$6-1}')  

          samtools faidx ${fasta} "$irb" |sed 1d > IRb.fa
          samtools faidx ${fasta} "$ira" |sed 1d > IRa.fa 
          samtools faidx ${fasta} "$ssc" |sed 1d > ssc.fa 
          samtools faidx ${fasta} "$lsc" |sed 1d > lsc.fa 

          cat header lsc.fa IRa.fa ssc.fa IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'| sed "s/^>.*$/>${prefix}/" >${out}
          rm _tmp.blast.out lsc.fa ssc.fa IRa.fa IRb.fa header ${fasta}.fai
          echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
        fi 
      fi  
    else 
          ## otherwise check if it is split at LSC or SSC
          end1=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print '$seq_len'-$3}')
          end2=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $6}')
          sum=$(($end1 + $end2))
          mid=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $2-$5}')

          if [ $mid -gt $sum ];then 
            ## the fasta is split at SSC region
            lsc=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$5+1"-"$2-1}')
            ira=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$2"-"$3}')
            ssc_start=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$3+1"-"'$seq_len'}')
            ssc_end=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":1-"$6-1}')
            irb=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$6"-"$5}')

            samtools faidx ${fasta} "$irb" |sed 1d > IRb.fa
            samtools faidx ${fasta} "$ira" |sed 1d > IRa.fa 
            samtools faidx ${fasta} "$ssc_start" |sed 1d > sscs.fa 
            samtools faidx ${fasta} "$ssc_end" |sed 1d > ssce.fa 
            samtools faidx ${fasta} "$lsc" |sed 1d > lsc.fa 

            cat header lsc.fa IRa.fa sscs.fa ssce.fa IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'| sed "s/^>.*$/>${prefix}/" > ${out} 
            rm _tmp.blast.out lsc.fa sscs.fa ssce.fa IRa.fa IRb.fa header ${fasta}.fai
            echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"

          else 
                ## the fasta is split at LSC region
                lsc_start=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$3+1"-"'$seq_len'}')
                lsc_end=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":1-"$6-1}')
                ira=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$6"-"$5}')
                ssc=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$5+1"-"$2-1}')
                irb=$(cat _tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$2"-"$3}')

                samtools faidx ${fasta} "$irb" |sed 1d > IRb.fa
                samtools faidx ${fasta} "$ira" |sed 1d > IRa.fa 
                samtools faidx ${fasta} "$lsc_start" |sed 1d > lscs.fa 
                samtools faidx ${fasta} "$lsc_end" |sed 1d > lsce.fa 
                samtools faidx ${fasta} "$ssc" |sed 1d > ssc.fa 

                cat header lscs.fa lsce.fa IRa.fa ssc.fa IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'| sed "s/^>.*$/>${prefix}/" > ${out} 
                rm _tmp.blast.out lscs.fa lsce.fa ssc.fa IRa.fa IRb.fa header ${fasta}.fai
                echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
          fi
    fi
fi

cd ../
rm -r _tmp_$(basename $fasta)
