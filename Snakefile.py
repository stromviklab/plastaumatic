import os
import glob

CULTIVAR,DIRECTION= glob_wildcards("/home/norac97/scratch/Glycine/{cultivar}_{direction}.fastq.gz")
OUTCOME= ["1","2"]

start1insert="{print $4}"
end1insert="{print $5}"
end2insert="{print $5}"

rule all:
    input:
        expand("trimmedReads/{cultivar}_{direction}P.fastq", cultivar=CULTIVAR, direction=DIRECTION),
        expand("assembledGenome/QUAST_OUT_{outcome}_{cultivar}/report.txt", cultivar=CULTIVAR, outcome=OUTCOME),
        expand("annotatedGenome/{outcome}_{cultivar}_full.gb", cultivar=CULTIVAR, outcome=OUTCOME)
        
rule decompress:
    input:
        compressed="/home/norac97/scratch/Glycine/{cultivar}_{direction}.fastq.gz"
    output:
        decompressed="/home/norac97/scratch/Glycine/{cultivar}_{direction}.fastq"
    threads:
        5
    shell:
        """
        pigz -d -p 5 -k {input.compressed}
        """

rule trimmomatic:
    input:
        read1="/home/norac97/scratch/Glycine/{cultivar}_R1.fastq",
        read2="/home/norac97/scratch/Glycine/{cultivar}_R2.fastq"
    output:
        forwardPaired="trimmedReads/{cultivar}_R1P.fastq",
        reversePaired="trimmedReads/{cultivar}_R2P.fastq", 
        forwardUnPaired="trimmedReads/{cultivar}_R1U.fastq", 
        reverseUnPaired="trimmedReads/{cultivar}_R2U.fastq"
    threads:
        10
    params:
        basename="trimmedReads/{cultivar}.fastq",
        log="trimmedReads/{cultivar}.log"
    shell:
        """
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads {threads} {input.read1} {input.read2} {output.forwardPaired} {output.forwardUnPaired} {output.reversePaired} {output.reverseUnPaired} ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60
        """

rule novoplasty:
    input:
        configfile="config_{cultivar}.txt",
        fowardPaired="trimmedReads/{cultivar}_R1P.fastq",
        reversePaired="trimmedReads/{cultivar}_R2P.fastq"
    output:
        option1="assembledGenome/Option_1_{cultivar}.fasta",
        option2="assembledGenome/Option_2_{cultivar}.fasta"
    shell:
        """
        perl /home/norac97/scratch/NOVOPlasty/NOVOPlasty4.3.1.pl -c {input.configfile}
        """

rule quast:
    input:
        option="assembledGenome/Option_{outcome}_{cultivar}.fasta"
    output:
        quastout="assembledGenome/QUAST_OUT_{outcome}_{cultivar}/report.txt"
    threads:
        5
    shell:
        """
        quast.py -t 5 -o assembledGenome/QUAST_OUT_{wildcards.outcome}_{wildcards.cultivar} {input.option}
        """

rule nucmer:
    input:
        ref="Glycine_max.fasta",
        option="assembledGenome/Option_{outcome}_{cultivar}.fasta"
    output:
        delta="prep/delta_{outcome}_{cultivar}.delta",
        coords="prep/coord_{outcome}_{cultivar}.coords"
    threads:
        5
    shell:
        """
        nucmer -t 5 --delta={output.delta} {input.ref} {input.option}
        """
        """
        show-coords {output.delta} > {output.coords}
        """

rule samtools:
    input:
        option="assembledGenome/Option_{outcome}_{cultivar}.fasta",
        coords="prep/coord_{outcome}_{cultivar}.coords"
    output:
        region1="prep/{outcome}_{cultivar}_region1.fa",
        region2="prep/{outcome}_{cultivar}_region2.fa",
        full="standardizedGenome/{outcome}_{cultivar}_full.fa"
    shell:
        """
        samtools faidx {input.option}
        """
        """
        samtools faidx {input.option} Contig1:$(cat {input.coords} |sed '1,5d' |awk '$1==1 {start1insert}')-$(cat {input.coords} |sed '1,5d' |awk '$5>150000 {end1insert}') > {output.region1}
        """
        """
        samtools faidx {input.option} Contig1:0-$(cat {input.coords} |sed '1,5d' |awk '$4<=1 {end2insert}') |sed 1d > {output.region2}
        """
        """
        cat {output.region1} {output.region2} > {output.full}
        """

rule PGA:
    input:
        full="standardizedGenome/{outcome}_{cultivar}_full.fa"
    output:
        outfiles="annotatedGenome/{outcome}_{cultivar}_full.gb"
    shell:
        """
        perl /home/norac97/scratch/PGA/PGA.pl -r GBref -t standardizedGenome -o annotatedGenome
        """
