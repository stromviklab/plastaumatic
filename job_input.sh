#!/bin/bash

dir=/home/norac97/scratch/Glycine #Input here the path to the directory containing raw reads
current_dir=/home/norac97/scratch/automated_pipeline_norac97 
config_base=/home/norac97/scratch/automated_pipeline_norac97/config_novo.txt 
snakefile_base=/home/norac97/scratch/automated_pipeline_norac97/Snakefile_single.py
seed=/home/norac97/scratch/automated_pipeline_norac97/G_max_seed_atpA.fasta
ref_fasta=/home/norac97/scratch/automated_pipeline_norac97/Glycine_max.fasta
ref_gb=/home/norac97/scratch/automated_pipeline_norac97/GBref
max_memory=40
threads_available=16
path_to_trimmomatic=$EBROOTTRIMMOMATIC/trimmomatic-0.39.jar
path_to_novoplasty=/home/norac97/scratch/NOVOPlasty/NOVOPlasty4.3.1.pl
path_to_PGA=/home/norac97/scratch/PGA/PGA.pl

mkdir all_snakefiles

for entry in "${dir}"/*
do
  mv $entry $(echo "$entry" |sed -r 's|_|-|g; s|-([^-]*)$|\_\1|')
done

for entry in "${dir}"/*
do
  mkdir $(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)
  cat ${config_base} |sed "s|test|$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)|g" |sed "s|max_memory|${max_memory}|g" |sed "s|path_to_seed|${seed}|g" |sed "s|path_to_reference|${ref_fasta}|g" |sed "s|each_dir|${current_dir}/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)|g" > $(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_config.txt
done

for entry in "${dir}"/*
do
  if [[ "$(echo "$entry" | rev | cut -d '_' -f 1 | rev)" == *1* ]]; then
    cat $(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_config.txt |sed "s|read1|$(basename -s .fastq "${entry%.*}")P.fastq|g" > $(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_config.tmp && mv $(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_config.tmp $(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_config.txt
    cat ${snakefile_base} |sed "s|forward_read|$(basename -s .fastq "${entry%.*}")|g" > all_snakefiles/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_snakefile.py
  fi
done

for entry in "${dir}"/*
do
  if [[ "$(echo "$entry" | rev | cut -d '_' -f 1 | rev)" == *2* ]]; then
    cat $(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_config.txt |sed "s|read2|$(basename -s .fastq "${entry%.*}")P.fastq|g" > $(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_config.tmp && mv $(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_config.tmp $(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_config.txt
    cat all_snakefiles/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_snakefile.py |sed "s|reverse_read|$(basename -s .fastq "${entry%.*}")|g" > all_snakefiles/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_snakefile.tmp && mv all_snakefiles/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_snakefile.tmp all_snakefiles/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_snakefile.py
    cat all_snakefiles/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_snakefile.py |sed "s|path_to_ref_fasta|${ref_fasta}|g" |sed "s|path_to_ref_gb|${ref_gb}|g" |sed "s|path_to_trimmomatic|${path_to_trimmomatic}|g" |sed "s|path_to_novoplasty|${path_to_novoplasty}|g" |sed "s|path_to_PGA|${path_to_PGA}|g" |sed "s|dir_to_raw_reads|${dir}|g" |sed "s|dir_to_each_assembly|${current_dir}/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)|g" |sed "s|reads_basename|$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)|g" > all_snakefiles/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_snakefile.tmp && mv all_snakefiles/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_snakefile.tmp all_snakefiles/$(basename -s .fastq "${entry%.*}" | cut -d '_' -f 1)_snakefile.py
  fi
done

for entry in "all_snakefiles"/*
do
  cat $current_dir/run.sh |sed "s|each_snakefile|all_snakefiles/$(basename "$entry")|g" > $(basename "${entry%.*}")_job.sh
  sbatch $(basename "${entry%.*}")_job.sh
done