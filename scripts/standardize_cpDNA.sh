#!/bin/bash
set -eo pipefail
while getopts 'd:g:o:p:' options
do
  case $options in
    d) dir=$OPTARG ;;
    g) gb=$OPTARG ;;
    o) out=$OPTARG ;;
    p) prefix=$OPTARG ;;
  esac
done


standardize() {
  echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tstandardizing the input fasta file"

  fasta=$1
  option=$2
  ## index the fasta input and get the sequence length
  samtools faidx ${fasta}
  seq_len=$(cut -f2 ${fasta}.fai)

  mkdir -p ${dir}/_tmp_plastaumatic 
  dir2=${dir}/_tmp_plastaumatic

  ## get the header 
  head -1 ${fasta} > ${dir2}/header

  ## blast the fasta to itself to get the repeat alignments
  ## ideally, the IR sequences in plastome should be identical, and hence a perc_identity of 99 is used
  ## self match is removed and any match less than 1000 bp length is removed
  ## based on test runs, the length can be reduced down to 100 bp (not sure how it preforms with 100 though)
  blastn -query ${fasta} -subject ${fasta} -perc_identity 99 -evalue 0.00001 -outfmt '6 qseqid qstart qend sseqid sstart send length pident'|awk '$2!=$5&&$3!=$6 {print}'|sort  -k7,7nr -k2,2n  > ${dir2}/_tmp.blast.out

  head -2 ${dir2}/_tmp.blast.out > ${dir2}/_tmp.blast.out2 
  start=$(head -1 ${dir2}/_tmp.blast.out|awk '{print $2}')
  end=$(head -1 ${dir2}/_tmp.blast.out|awk '{print $5}')
  if [ $start -eq 1 ];then  
    awk '$3=='$seq_len' && $6-1=='$end' {print}' ${dir2}/_tmp.blast.out >> ${dir2}/_tmp.blast.out2 
    awk '$5=='$seq_len' && $2-1=='$end' {print}' ${dir2}/_tmp.blast.out >> ${dir2}/_tmp.blast.out2 
  elif [ $end -eq $seq_len ];then
    awk '$2==1 && $5+1=='$start' {print}' ${dir2}/_tmp.blast.out >> ${dir2}/_tmp.blast.out2
    awk '$6==1 && $3+1=='$start' {print}' ${dir2}/_tmp.blast.out >> ${dir2}/_tmp.blast.out2
  fi
  mv ${dir2}/_tmp.blast.out2 ${dir2}/_tmp.blast.out

  ## if the fasta is split at inverted repeats, it would give 4 alignments
  if [ $(grep -c "^" ${dir2}/_tmp.blast.out) -gt 2 ];then 
      ir_tmp=$(cat ${dir2}/_tmp.blast.out|awk '$3=='$seq_len' {print $6-1}')
      sc1=$(cat ${dir2}/_tmp.blast.out|awk '$5=='$ir_tmp' {print $6-$3}')
      sc2=$(cat ${dir2}/_tmp.blast.out|awk '$3=='$seq_len' {print $2-$5}')
      if [ $sc2 -gt $sc1 ];then 
          ## do this if the fasta is split at IRa region 
          lsc=$(cat ${dir2}/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$5+1"-"$2-1}')
          ira_start=$(cat ${dir2}/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$2"-"$3}')
          ira_end=$(cat ${dir2}/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$2"-"$3}') 
          ssc=$(cat ${dir2}/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$3+1"-"$6-1}')
          irb_start=$(cat ${dir2}/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$6"-"$5}') 
          irb_end=$(cat ${dir2}/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$6"-"$5}')
      else 
          ## do this if the fasta is split at IRb region 
          lsc=$(cat ${dir2}/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$3+1"-"$6-1}')
          ira_start=$(cat ${dir2}/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$6"-"$5}') 
          ira_end=$(cat ${dir2}/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$6"-"$5}')
          ssc=$(cat ${dir2}/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$5+1"-"$2-1}')
          irb_start=$(cat ${dir2}/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$2"-"$3}')
          irb_end=$(cat ${dir2}/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$2"-"$3}') 
      fi 
      samtools faidx ${fasta} "$irb_start" |sed "s/$irb_start/IRb/" > ${dir2}/irbs 
      samtools faidx ${fasta} "$irb_end" |sed 1d > ${dir2}/irbe  
      samtools faidx ${fasta} "$ira_start" |sed "s/$ira_start/IRa/" > ${dir2}/iras 
      samtools faidx ${fasta} "$ira_end" |sed 1d > ${dir2}/irae  
      cat ${dir2}/irbs ${dir2}/irbe > ${dir2}/IRb.fa 
      cat ${dir2}/iras ${dir2}/irae > ${dir2}/IRa.fa 
      ## double checking the IR alignments 
      if [ "$(blastn -query ${dir2}/IRa.fa -subject ${dir2}/IRb.fa -perc_identity 99 -evalue 0.00001 -outfmt '6 qseqid qstart qend sseqid sstart send length pident' |awk 'NR==1&&$2==$6&&$3==$5 {print "Success"}')" == "Success" ] ;then 
        samtools faidx ${fasta} "$lsc" |sed 1d > ${dir2}/LSC.fa 
        samtools faidx ${fasta} "$ssc" |sed 1d > ${dir2}/SSC.fa 
        cat ${dir2}/IRa.fa | sed 1d > ${dir2}/_tmp.IRa.fa 
        cat ${dir2}/IRb.fa | sed 1d > ${dir2}/_tmp.IRb.fa
        cat ${dir2}/header ${dir2}/LSC.fa ${dir2}/_tmp.IRa.fa ${dir2}/SSC.fa ${dir2}/_tmp.IRb.fa | sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' | sed "s/^>.*$/>${prefix}/" > ${dir2}/out_${option}.fa
        echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
      fi

  else 
      ## if the end or start the sequence is there in the alignemnts, that means the IR locations needs to be adjusted
      if [ "$(cat ${dir2}/_tmp.blast.out| awk '$3=='$seq_len'||$2==1 {print "yes"}')" == "yes" ];then 
        if [ "$(cat ${dir2}/_tmp.blast.out| awk '$3=='$seq_len' {print "yes"}')" == "yes" ];then
          if [ $(awk '$3=='$seq_len' {print $6}' ${dir2}/_tmp.blast.out) -gt 50000 ];then
            sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' ${fasta} | sed "s/^>.*$/>${prefix}/" > ${dir2}/out_${option}.fa
            echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
          else 
            ira=$(cat ${dir2}/_tmp.blast.out |awk '$3=='$seq_len' {print $1":"$2"-"$3}')
            irb=$(cat ${dir2}/_tmp.blast.out |awk '$3=='$seq_len' {print $1":"$6"-"$5}')
            lsc=$(cat ${dir2}/_tmp.blast.out |awk '$3=='$seq_len' {print $1":"$5+1"-"$2-1}')  
            ssc=$(cat ${dir2}/_tmp.blast.out |awk '$3=='$seq_len' {print $1":1-"$6-1}')  

            samtools faidx ${fasta} "$irb" |sed 1d > ${dir2}/IRb.fa
            samtools faidx ${fasta} "$ira" |sed 1d > ${dir2}/IRa.fa 
            samtools faidx ${fasta} "$ssc" |sed 1d > ${dir2}/ssc.fa 
            samtools faidx ${fasta} "$lsc" |sed 1d > ${dir2}/lsc.fa 

            cat ${dir2}/header ${dir2}/lsc.fa ${dir2}/IRa.fa ${dir2}/ssc.fa ${dir2}/IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' | sed "s/^>.*$/>${prefix}/" > ${dir2}/out_${option}.fa
            echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
          fi 
        fi 
        if [ "$(cat ${dir2}/_tmp.blast.out| awk '$2==1 {print "yes"}')" == "yes" ];then
          if [ $(cat ${dir2}/_tmp.blast.out |awk '$2==1 {print '$seq_len'-$5}') -gt 50000 ];then 
            ira=$(cat ${dir2}/_tmp.blast.out |awk '$2==1 {print $1":"$2"-"$3}')
            irb=$(cat ${dir2}/_tmp.blast.out |awk '$2==1 {print $1":"$6"-"$5}')
            lsc=$(cat ${dir2}/_tmp.blast.out |awk '$2==1 {print $1":"$5+1"-"'$seq_len'}')  
            ssc=$(cat ${dir2}/_tmp.blast.out |awk '$2==1 {print $1":"$3+1"-"$6-1}')  

            samtools faidx ${fasta} "$irb" |sed 1d > ${dir2}/IRb.fa
            samtools faidx ${fasta} "$ira" |sed 1d > ${dir2}/IRa.fa 
            samtools faidx ${fasta} "$ssc" |sed 1d > ${dir2}/ssc.fa 
            samtools faidx ${fasta} "$lsc" |sed 1d > ${dir2}/lsc.fa 

            cat ${dir2}/header ${dir2}/lsc.fa ${dir2}/IRa.fa ${dir2}/ssc.fa ${dir2}/IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'| sed "s/^>.*$/>${prefix}/" > ${dir2}/out_${option}.fa
            echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
          else 
            irb=$(cat ${dir2}/_tmp.blast.out |awk '$2==1 {print $1":"$2"-"$3}')
            ira=$(cat ${dir2}/_tmp.blast.out |awk '$2==1 {print $1":"$6"-"$5}')
            ssc=$(cat ${dir2}/_tmp.blast.out |awk '$2==1 {print $1":"$5+1"-"'$seq_len'}')  
            lsc=$(cat ${dir2}/_tmp.blast.out |awk '$2==1 {print $1":"$3+1"-"$6-1}')  

            samtools faidx ${fasta} "$irb" |sed 1d > ${dir2}/IRb.fa
            samtools faidx ${fasta} "$ira" |sed 1d > ${dir2}/IRa.fa 
            samtools faidx ${fasta} "$ssc" |sed 1d > ${dir2}/ssc.fa 
            samtools faidx ${fasta} "$lsc" |sed 1d > ${dir2}/lsc.fa 

            cat ${dir2}/header ${dir2}/lsc.fa ${dir2}/IRa.fa ${dir2}/ssc.fa ${dir2}/IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'| sed "s/^>.*$/>${prefix}/" > ${dir2}/out_${option}.fa
            echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
          fi 
        fi  
      else 
            ## otherwise check if it is split at LSC or SSC
            end1=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print '$seq_len'-$3}')
            end2=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $6}')
            sum=$(($end1 + $end2))
            mid=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $2-$5}')

            if [ $mid -gt $sum ];then 
              ## the fasta is split at SSC region
              lsc=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$5+1"-"$2-1}')
              ira=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$2"-"$3}')
              ssc_start=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$3+1"-"'$seq_len'}')
              ssc_end=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":1-"$6-1}')
              irb=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$6"-"$5}')

              samtools faidx ${fasta} "$irb" |sed 1d > ${dir2}/IRb.fa
              samtools faidx ${fasta} "$ira" |sed 1d > ${dir2}/IRa.fa 
              samtools faidx ${fasta} "$ssc_start" |sed 1d > ${dir2}/sscs.fa 
              samtools faidx ${fasta} "$ssc_end" |sed 1d > ${dir2}/ssce.fa 
              samtools faidx ${fasta} "$lsc" |sed 1d > ${dir2}/lsc.fa 

              cat ${dir2}/header ${dir2}/lsc.fa ${dir2}/IRa.fa ${dir2}/sscs.fa ${dir2}/ssce.fa ${dir2}/IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'| sed "s/^>.*$/>${prefix}/" > ${dir2}/out_${option}.fa
              echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"

            else 
                  ## the fasta is split at LSC region
                  lsc_start=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$3+1"-"'$seq_len'}')
                  lsc_end=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":1-"$6-1}')
                  ira=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$6"-"$5}')
                  ssc=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$5+1"-"$2-1}')
                  irb=$(cat ${dir2}/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$2"-"$3}')

                  samtools faidx ${fasta} "$irb" |sed 1d > ${dir2}/IRb.fa
                  samtools faidx ${fasta} "$ira" |sed 1d > ${dir2}/IRa.fa 
                  samtools faidx ${fasta} "$lsc_start" |sed 1d > ${dir2}/lscs.fa 
                  samtools faidx ${fasta} "$lsc_end" |sed 1d > ${dir2}/lsce.fa 
                  samtools faidx ${fasta} "$ssc" |sed 1d > ${dir2}/ssc.fa 

                  cat ${dir2}/header ${dir2}/lscs.fa ${dir2}/lsce.fa ${dir2}/IRa.fa ${dir2}/ssc.fa ${dir2}/IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'| sed "s/^>.*$/>${prefix}/" > ${dir2}/out_${option}.fa
                  echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tFinished"
            fi
      fi
  fi
}

if [ -f ${dir}/Circularized_assembly_1_${prefix}.fasta ];then
  standardize "${dir}/Circularized_assembly_1_${prefix}.fasta" 0
  mv ${dir2}/out_0.fa ${out}
  rm -r ${dir2}
elif [[ -f ${dir}/Option_1_${prefix}.fasta ]] && [[ -f ${dir}/Option_2_${prefix}.fasta ]] && [[ ! -f ${dir}/Option_3_${prefix}.fasta ]];then 
  opt1_len=$(cat ${dir}/Option_1_${prefix}.fasta |sed 1d|tr -d '\n'|wc -c)
  opt2_len=$(cat ${dir}/Option_2_${prefix}.fasta |sed 1d|tr -d '\n'|wc -c)
  if [ $opt1_len -eq $opt2_len ];then 
    standardize "${dir}/Option_1_${prefix}.fasta" 1
    standardize "${dir}/Option_2_${prefix}.fasta" 2
    # get reference fasta
    echo ">ref" > ${dir2}/ref.fa
    awk '/ORIGIN/,/\/\//' ${gb} |grep -v "ORIGIN"|grep -v "\//"|tr -d "[0-9]"|tr -d ' '|tr -d '\n' >> ${dir2}/ref.fa 
    ## Selecting option1 or option2 based on the reference
    for k in 1 2; do
      blastn -query ${dir2}/ref.fa -subject ${dir2}/out_${k}.fa -perc_identity 95 -evalue 0.1 -outfmt '6 qseqid qstart qend sseqid sstart send length sstrand qcovhsp'|awk '$7>1000 {print}' |grep "plus" > ${dir2}/_tmp.blast2.out
        if [ $(awk '{sum+=$9;}END{print sum}' ${dir2}/_tmp.blast2.out) -gt 95 ];then
          mv ${dir2}/out_${k}.fa ${out}
        fi
    done
    if [ ! -f ${out} ];then 
      mv ${dir2}/out_1.fa ${out}
    fi 
    rm -r ${dir2}
  else 
    echo "ERROR: Two options from NOVOPlasty are of different size"
    exit 1
  fi 
elif [ -f ${dir}/Option_1_${prefix}.fasta ] && [[ ! -f ${dir}/Option_3_${prefix}.fasta ]];then 
  standardize "${dir}/Option_1_${prefix}.fasta" 1
  mv ${dir2}/out_1.fa ${out}
  rm -r ${dir2}
elif [ -f ${dir}/Option_2_${prefix}.fasta ] && [[ ! -f ${dir}/Option_3_${prefix}.fasta ]];then 
  standardize "${dir}/Option_2_${prefix}.fasta" 2
  mv ${dir2}/out_2.fa ${out}
  rm -r ${dir2}
else 
  echo "ERROR: NOVOPlasty did not return a complete circular assembly or too many assemblies"
  exit 1
fi 


## Check and correct non-ACTG characters in the assembly 
echo -e "$(date +'%Y-%m-%d %H:%M:%S')\tchecking for non-ATCG characters in the $prefix assembly"
non_actg=$(cat ${out}|sed 1d|tr -d '\n' |tr -d 'ACTGactg')
while [ ! -z "$non_actg" ];do 
  samtools faidx ${out}
  header=$(awk '{print $1}' ${out}.fai) 
  line1=$(cat ${out}|sed 1d|tr -d '\n' |sed 's/./&\n/g'|grep -n "[A-Za-z]"|grep -Ev "[ACTGactg]"|head -1|sed 's/:/\t/g')
  char=$(echo $line1|awk '{print $2}')
  if [ $(echo $line1|awk '{print $1}') -gt 100 ];then 
    loc=$(echo $line1|awk '{print "'$header':"$1-100"-"$1-1}') 
    seq=$(samtools faidx -n 200 ${out} "$loc"|sed 1d)
    new_seq=$(grep -o "$seq[A-Za-z]" ${dir}/Assembled_reads_${prefix}_R*.fasta |cut -d ':' -f2|sort|uniq -c |sort -k1 -n -r |head -1|awk '{print $2}')
    sed -i "s/$seq$char/$new_seq/1" ${out} 
    non_actg=$(cat ${out}|sed 1d|tr -d '\n' |tr -d 'ACTGactg')
  else
    loc=$(echo $line1|awk '{print "'$header':"$1+1"-"$1+100}')
    seq=$(samtools faidx -n 200 ${out} "$loc"|sed 1d)
    new_seq=$(grep -o "[A-Za-z]$seq" ${dir}/Assembled_reads_${prefix}_R*.fasta |cut -d ':' -f2|sort|uniq -c |sort -k1 -n -r |head -1|awk '{print $2}')
    sed -i "s/$char$seq/$new_seq/1" ${out} 
    non_actg=$(cat ${out}|sed 1d|tr -d '\n' |tr -d 'ACTGactg')
  fi 
done 
echo -e "$(date +'%Y-%m-%d %H:%M:%S')\t$prefix Plastome assembly size: $(cat ${out} |sed 1d|tr -d '\n'|wc -c) bps"
