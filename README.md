# Plastaumatic: An automated pipeline to assemble and annotate plastome sequences. 
![](./pipeline.pdf)

## Introduction
**`plastaumatic`** is an automated pipeline developed for both assembly and annotation of plastomes, with the scope of the researcher being able to load whole genome sequence data with minimal manual input, and therefore a faster runtime. The main structure of the current automated pipeline includes trimming of adaptor and low-quality sequences using `fastp`, *de novo* plastome assembly using `NOVOPlasty`, standardization and quality checking of the assembled genomes through a custom script utilizing `BLAST+` and `SAMtools`, annotation of the assembled genomes using `AnnoPlast`, and finally generating required files for NCBI GenBank submissions.

This repository includes a `snakefile` which uses `Snakemake` workflow manager to perform all the tasks. Also, a shell executable that does the same. 

### Pre-requisites:
[**`Snakemake v5.10.0+`**](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html) [if using the `Snakemake` pipeline] </br>
[**`fastp v0.23.0`**](https://github.com/OpenGene/fastp) </br>
[**`NOVOPlasty v4.3.1`**](https://github.com/ndierckx/NOVOPlasty)  </br>
**`Python v3.8`**  </br>
**`BioPython v1.79` & `pandas v1.4.3`** (for [**`AnnoPlast`**](./scripts/AnnoPlast.py)) </br>
[**`Blast+ v2.12.0`**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) </br>
[**`samtools v1.9`**](http://www.htslib.org/download/) </br>


### installation
Simply clone this repository in your desired destination and add to your PATH</br>

    git clone https://github.com/stromviklab/plastaumatic.git
    export PATH=$(pwd)/plastaumatic:$PATH
    echo 'export PATH='$(pwd)'/plastaumatic:$PATH' >> ~/.bashrc


### how-to-run (Snakemake)
1. Before running the Snakemake pipeline, make sure all the pre-requisites are installed and available in your PATH

2. Make sure to edit the `config.yaml` file before your run (a template `config.yaml` file is available from the repository)

````
seed: Test_dataset/seed.fasta                                       # path to a seed fasta file for NOVOPlasty (see: https://github.com/ndierckx/NOVOPlasty)
ref_gb: Test_dataset/reference.gb                                   # path to a reference GenBank file for annotation (see: https://github.com/SaiReddy-A/AnnoPlast) 
range: 140000-160000                                                # estimated plastome size range 
repo: plastaumatic                                                  # path to the plastaumatic repository
novo_path: /mnt2/software/NOVOPlasty/NOVOPlasty4.3.1.pl             # path to the NOVOPlasty executable
samples:                                                            # prefix and reads of one/many genomes go under samples
  genome1: Test_dataset/test1.fq.gz,Test_dataset/test2.fq.gz            # prefix: forward_read,revese_read
````
* If you want to run this pipeline on multiple genomes, just add more lines below `samples` in this format:
````
samples:
    genome1: forward.fq,reverse.fq
    genome2: forward.fq,reverse.fq
    genome3: forward.fq,reverse.fq
    ...
    ...
    genomeX: forward.fq,reverse.fq
````

3. Copy the `snakefile` to your preferred output directory. Modify the path to config file in the `snakefile` and run </br>

````
snakemake
````    


### how-to-run (shell script)
1. Before running the `plastaumatic` executable, make sure all the pre-requisites are installed and available in your PATH
2. run `plastaumatic -h` for help</br>
````
Usage: plastaumatic -s seed.fa -g reference.gb -r <140000-160000> -f fof.txt -n NOVOPlasty4.3.1.pl

options:
         -s      Path to the seed file for assembly</br>
         -g      Path to the reference GenBank file</br>
         -r      Plastome assembly size range [140000-160000]</br>
         -f      Path to the file-of-filenames with reads</br>
         -n      Path to NOVOPlasty executable</br>
         -h      Shows this help message</br>
````

* The file-of-filenames is a simple txt file with comma seperated prefix and read files of one/many genomes   
````
genome1,forward.fq,reverse.fq
genome2,forward.fq,reverse.fq
genome3,forward.fq,reverse.fq
````

3. An example `plastaumatic` run </br>

````
plastaumatic -s Test_dataset/seed.fasta -g Test_dataset/reference.gb -r 140000-160000 -f readList.txt -n software/NOVOPlasty/NOVOPlasty4.3.1.pl      
````    
> $ cat readList.txt</br>
> genome1,Test_dataset/test1.fq.gz,Test_dataset/test2.fq.gz




### Output files
For each genome a *prefix* directory is created from where the pipeline is run. Each *prefix* directory contains
```
prefix
├── 00-logs
├── 01-trim
├── 02-assemble
├── 03-standardize
├── 04-annotate
├── 05-tbl
├── prefix_config.txt
├── prefix.plastome.fa # a symlink to the final plastome assembly file  
└── prefix.plastome.gb # a symlink to the final plastome annotation file
```




### Citations
Since `plastaumatic` uses multiple software in its pipeline, publishing the results obtained form `plastaumtic` should also cite the following sources

> 1. Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
> 2. Dierckxsens, N., Mardulyn, P., & Smits, G. (2017). NOVOPlasty: de novo assembly of organelle genomes from whole genome data. Nucleic acids research, 45(4), e18. https://doi.org/10.1093/nar/gkw955
> 3. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC Bioinformatics. 2009 Dec 15;10:421. doi: 10.1186/1471-2105-10-421. PMID: 20003500; PMCID: PMC2803857
> 4. Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352


