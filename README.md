# Plastaumatic: An automated pipeline to assemble and annotate plastome sequences. 


## Introduction
**`plastaumatic`** is an automated pipeline developed for both assembly and annotation of plastomes, with the scope of the researcher being able to load whole genome sequence data with minimal manual input, and therefore a faster runtime. The main structure of the current automated pipeline includes file decompression by `pigz`, trimming of adaptor and low-quality sequences using `Trimmomatic`, *de novo* plastome assembly using `NOVOPlasty`, standardization of the assembled genomes through a custom script utilizing `BLAST+` and `SAMtools`, annotation of the assembled genomes using `PGA`, and a final inspection of the annotation using a custom script.

This pipeline uses `Snakemake` workflow manager to perform the analyses. However, a shell script is written to modify the input parameters in various steps and this is the main executable for running this pipeline. 

### Pre-requisites:
[**`Snakemake`**](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html) </br>
**`pigz`** </br>
[**`Trimmomatic`**](https://github.com/usadellab/Trimmomatic) </br>
[**`NOVOPlasty`**](https://github.com/ndierckx/NOVOPlasty) </br>
[**`PGA`**](https://github.com/quxiaojian/PGA) </br>
[**`Blast+`**](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) </br>
[**`samtools`**](http://www.htslib.org/download/) </br>


### installation
Simply clone this repository in your desired destination</br>

    git clone https://github.com/stromviklab/plastaumatic.git


### how-to-run
1. Before running the executable, make sure all the pre-requisites are installed properly

2. Make sure to edit the config file for the executable before your run (a template for config file is `test.cfg` in the installation directory)

````prefix=test                                                           #prefix to be used for the output directories and file names 
read1=/home/noracwy/pipeline-test/read1.fq                                #complete path to the forward reads (can be .fq or .fq.gz)
read2=/home/noracwy/pipeline-test/read2.fq                                #complete path to the reverse reads (can be .fq or .fq.gz)
seed=/home/noracwy/pipeline-test/seed.fasta                               #complete path to a seed fasta file for NOVOPlasty (see: https://github.com/ndierckx/NOVOPlasty)
ref_fasta=/home/noracwy/pipeline-test/ref.fasta                           #complete path to a reference fasta file for NOVOPlasty (see: https://github.com/ndierckx/NOVOPlasty)
ref_gb=/home/noracwy/pipeline-test/GBref/                                 #complete path to a directory with reference GenBank file for PGA (see: https://github.com/quxiaojian/PGA)
max_memory=40                                                             #maximum memory available to use in GB 
threads_available=1                                                       #maximum numer of threads available to use 
range=120000-160000                                                       #estimated plastome size range 
out=/path/to/outDir/                                                      #complete path to the output directory (a new directory named prefix will be created here)
path_to_trimmomatic=/home/noracwy/Trimmomatic-0.39/trimmomatic-0.39.jar   #complete path to the Trimmomatic jar file
adapters=/home/noracwy/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa          #complete path to the adapters fasta file
path_to_novoplasty=/home/noracwy/NOVOPlasty/NOVOPlasty4.3.1.pl            #complete path to the NOVOPlasty executable
path_to_PGA=/home/noracwy/PGA/PGA.pl                                      #complete path to the PGA executable
path_to_repo=/home/noracwy/pipeline-test/plastaumatic/                    #complete path to the plastaumatic repository
````


3. The pipeline can be executed by simply running  </br>

````console
run_plastaumatic.sh input.cfg
````    


### Output files
```
prefix
├── reads
├── trimmedReads
├── assembledGenome
├── standardizedGenome
├── annotatedGenome
├── logs
└── snakefile
```

A main directory named *prefix* is created in the specified *outDir* from the config file, and inside this prefix directory you can find reads and trimmed reads, NOVOPlasty output files in *assembledGenome*, and the standardized assembly fasta in *standardizedGenome*, and the annotations in *annotatedGenome*. 




### Citations
Since `plastaumatic` uses multiple software in its pipeline, publishing the results obtained form `plastaumtic` should also cite the following sources

> 1. Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics (Oxford, England), 30(15), 2114–2120. https://doi.org/10.1093/bioinformatics/btu170
> 2. Dierckxsens, N., Mardulyn, P., & Smits, G. (2017). NOVOPlasty: de novo assembly of organelle genomes from whole genome data. Nucleic acids research, 45(4), e18. https://doi.org/10.1093/nar/gkw955
> 3. Qu X-J, Moore MJ, Li D-Z, Yi T-S. 2019. PGA: a software package for rapid, accurate, and flexible batch annotation of plastomes. Plant Methods 15:50.
> 4. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC Bioinformatics. 2009 Dec 15;10:421. doi: 10.1186/1471-2105-10-421. PMID: 20003500; PMCID: PMC2803857
> 5. Heng Li, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, 1000 Genome Project Data Processing Subgroup, The Sequence Alignment/Map format and SAMtools, Bioinformatics, Volume 25, Issue 16, 15 August 2009, Pages 2078–2079, https://doi.org/10.1093/bioinformatics/btp352


