#!/bin/bash

if [ $# -ne 2 ];
then
    echo "Usage: 
        check_internal_stops.sh input.gb input.fa"
    exit 1
fi

if [ -f $1 ] && [ -f $2 ]
then 

wget -q https://raw.githubusercontent.com/david-a-parry/translateDna/master/translateDna.pl

samtools faidx $2
grep -A1 "CDS  " $1 |grep -v "\--"|paste -d '\t' - - |awk '{print $3,$2}' OFS='\t'|awk -F '"' '{print $2,$0}' |awk '{print $1,$3}' OFS='\t' > _list

id=$(cut -f1 $2.fai)
echo -e "CHECKING FOR INTERNAL STOP CODONS ........."
grep -E -v "join|complement" _list |grep -v -wF "rps12"|while read i 
do
name=$(echo $i|awk '{print $1}')
loc=$(echo $i|awk '{print $2}'|awk -F '.' '{print "'$id':"$1"-"$3}')

count=$(samtools faidx  $2 $loc|grep -v ">" |perl translateDna.pl -i -|grep -c '*')
if [ $count -gt 1 ]
then 
ft=$(echo $i|awk '{print $2}')
echo -e "$name:$ft gene feature has internal stop codons"
fi 
done


grep "join" _list|grep -v "complement" |grep -v -wF "rps12"|while read i
do
    name=$(echo $i|awk '{print $1}')
    loc=$(echo "$i"|cut -f2|tr -d 'join()'|sed 's/\../-/g'|tr ',' '\n'|awk '{print "'$id':"$1}')
    count=$(samtools faidx $2 `echo $loc` |grep -v '>'|tr -d '\n'|perl translateDna.pl -i -|grep -c "*")
    
    if [ $count -gt 1 ];then 
    ft=$(echo $i|awk '{print $2}')
    echo -e "$name:$ft gene feature has internal stop codons"
    fi 
done

grep "complement" _list|grep -v "join"|grep -v -wF "rps12"|while read i 
do 
    name=$(echo $i|awk '{print $1}')
    loc=$(echo "$i"|cut -f2|tr -d 'complement()'|sed 's/\../-/g'|awk '{print "'$id':"$1}')
    count=$(samtools faidx -i $2 `echo $loc` |perl translateDna.pl -i -|grep -c "*")
    
    if [ $count -gt 1 ];then 
    ft=$(echo $i|awk '{print $2}')
    echo -e "$name:$ft gene feature has internal stop codons"
    fi 
done

grep "complement" _list|grep "join"|grep -v -wF "rps12"|while read i
do
    name=$(echo $i|awk '{print $1}')
    if [ $(echo $i|awk '{print $2}'|tr '(' '\t'|tr -d ')'|awk '{ if (NF==3&&$1=="complement"&&$2=="join") print "yes";else print "no"}') == "yes" ];then 

        loc=$(echo $i|awk '{print $2}'|tr '(' '\t'|tr -d ')'|awk 'NF==3&&$1=="complement"&&$2=="join" {print $3}'|sed 's/\../-/g'|tr ',' '\n'|sort -nr|awk '{print "'$id':"$1}')
        count=$(samtools faidx -i $2 `echo $loc` |grep -v ">"|tr -d '\n'|perl translateDna.pl -i -|grep -c "*")

        if [ $count -gt 1 ];then 
        ft=$(echo $i|awk '{print $2}')
        echo -e "$name:$ft gene feature has internal stop codons"
        fi 
    else 
    echo -e "Not able to determine the arrangement of $name CDS sequence"
    fi 

done

if [ $(grep -c -wF "rps12" _list) -eq 3 ];then 
    c1=$(grep -wF "rps12" _list |grep -v "join" |cut -f2 |tr -d 'complement()'|sed 's/\../-/g'|awk '{print "'$id':"$1}')
    c2=$(grep -wF "rps12" _list |grep "complement"|grep "join" |cut -f2 |tr -d 'complement(join)'|sed 's/\../-/g'|tr ',' '\n'|sed -n "2p"|awk '{print "'$id':"$1}')
    c3=$(grep -wF "rps12" _list |grep "complement"|grep "join" |cut -f2 |tr -d 'complement(join)'|sed 's/\../-/g'|tr ',' '\n'|sed -n "1p"|awk '{print "'$id':"$1}')
    count=$(samtools faidx -i $2 `echo $c1 $c2 $c3` |grep -v ">"|tr -d '\n'|perl translateDna.pl -i -|grep -c "*")
     
        if [ $count -gt 1 ];then 
        ft=$(grep -wF "rps12" _list |grep "complement"|grep "join" |awk '{print $2}')
        echo -e "rps12:$ft gene feature has internal stop codons"
        fi 

    j2=$(grep -wF "rps12" _list |grep -v "complement"|grep "join" |cut -f2 |tr -d 'join()'|sed 's/\../-/g'|tr ',' '\n'|sed -n "1p"|awk '{print "'$id':"$1}')
    j3=$(grep -wF "rps12" _list |grep -v "complement"|grep "join" |cut -f2 |tr -d 'join()'|sed 's/\../-/g'|tr ',' '\n'|sed -n "2p"|awk '{print "'$id':"$1}')
    j2_3=$(samtools faidx $2 `echo $j2 $j3` |grep -v ">"|tr -d '\n')
    c1fa=$(samtools faidx -i $2 `echo $c1` |grep -v ">"|tr -d '\n')
    count=$(echo $c1fa $j2_3 |tr -d ' ' |perl translateDna.pl -i -|grep -c "*")

        if [ $count -gt 1 ];then 
        ft=$(grep -wF "rps12" _list |grep -v "complement"|grep "join" |awk '{print $2}')
        echo -e "rps12:$ft gene feature has internal stop codons"
        fi 

else 
    echo -e "Not able to determine the arrangement of $name CDS sequences"

fi
rm translateDna.pl _list

else 
  echo "input files do not exist"

fi
