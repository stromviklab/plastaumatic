#!/bin/bash

if [ $# -ne 1 ];then
  echo "Usage: 
      loop.sh -d /path/to/dir/with-reads/"
  exit 1
fi

while getopts 'd:' options
do
  case $options in
    d) dir=$OPTARG ;;
  esac
done

if [ ! -d $dir ];then 
  echo "input directory does not exist"
  exit 1
fi 


if [ $(ls ${dir}/*fastq.gz|wc -l) > 0 ];then
    file_ext=fastq.gz
elif [ $(ls ${dir}/*fq.gz|wc -l) > 0 ];then
    file_ext=fq.gz
elif [ $(ls ${dir}/*.fastq|wc -l) > 0 ];then
    file_ext=fastq
elif [ $(ls ${dir}/*.fq|wc -l) > 0 ];then
    file_ext=fq
else 
    echo "ERROR: Can not determine the read file extension/ No read files exist"
fi 

declare -A prefix_array
base_name=$(ls ${dir}/*${file_ext}|awk -F '/' '{print $NF}'|sed "s/\.$file_ext//g")
if [ $base_name|awk -F '_' '{print $NF}' ]

for f in ${dir}/*${file_ext}; do
    prefix=basename $f .${file_ext}
    if [[ -z ${prefix_array[$prefix]} ]]; then 
            prefix_array[$prefix]=$f 
        else
            prefix_array[$prefix]+=" $f"
        fi
    done

    for prefix in "${!prefix_array[@]}"; do
        read -r f1 f2 <<< "${prefix_array[$prefix]}" 

        cat ${path_to_repo}/loop.cfg |sed "s|forward_file|"$f1"|g"| sed "s|reverse_file|"$f2"|" > "$prefix"_single.cfg
        chmod 777 "$prefix"_single.cfg
        bash ${path_to_repo}/run_plastaumatic.sh "$prefix"_single.cfg
    done

