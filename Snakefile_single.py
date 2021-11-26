wildcard_constraints:
    outcome="\d+"

start1insert="{print $4}"
end1insert="{print $5}"
end2insert="{print $5}"

OUTCOME= ["1","2"]

rule all:
    input:
        expand("dir_to_each_assembly/assembledGenome/QUAST_OUT_{outcome}_reads_basename/report.txt", outcome=OUTCOME),
        expand("dir_to_each_assembly/annotatedGenome/{outcome}_reads_basename_full.gb", outcome=OUTCOME)
        
rule decompress:
    input:
        compressed_forward="dir_to_raw_reads/forward_read.fastq.gz", 
        compressed_reverse="dir_to_raw_reads/reverse_read.fastq.gz"
    output:
        decompressed_forward="dir_to_raw_reads/forward_read.fastq",
        decompressed_reverse="dir_to_raw_reads/reverse_read.fastq"
    threads:
        threads_available
    shell:
        """
        echo "$(date): Starts decompressing reads_basename!"
        """
        """
        pigz -d -p {threads} -k {input.compressed_forward} {input.compressed_reverse}
        """
        """
        echo "$(date): Finished decompressing reads_basename!"
        """

rule trimmomatic:
    input:
        read1="dir_to_raw_reads/forward_read.fastq",
        read2="dir_to_raw_reads/reverse_read.fastq"
    output:
        forwardPaired="dir_to_each_assembly/trimmedReads/forward_readP.fastq",
        reversePaired="dir_to_each_assembly/trimmedReads/reverse_readP.fastq", 
        forwardUnPaired="dir_to_each_assembly/trimmedReads/forward_readU.fastq", 
        reverseUnPaired="dir_to_each_assembly/trimmedReads/reverse_readU.fastq"
    threads:
        threads_available
    shell:
        """
        echo "$(date): Trimmomatic for reads_basename is starting!"
        """
        """
        java -jar path_to_trimmomatic PE -threads {threads} {input.read1} {input.read2} {output.forwardPaired} {output.forwardUnPaired} {output.reversePaired} {output.reverseUnPaired} ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60
        """
        """
        echo "$(date): Trimmomatic for reads_basename is finished!"
        """

rule novoplasty:
    input:
        config_novo="dir_to_each_assembly/reads_basename_config.txt",
        fowardPaired="dir_to_each_assembly/trimmedReads/forward_readP.fastq",
        reversePaired="dir_to_each_assembly/trimmedReads/reverse_readP.fastq"
    output:
        option1="dir_to_each_assembly/assembledGenome/Option_1_reads_basename.fasta",
        option2="dir_to_each_assembly/assembledGenome/Option_2_reads_basename.fasta"
    shell:
        """
        echo "$(date): NOVOPlasty for reads_basename is starting!"
        """
        """
        perl path_to_novoplasty -c {input.config_novo}
        """
        """
        echo "$(date): NOVOPlasty for reads_basename is finished!"
        """

rule quast:
    input:
        option="dir_to_each_assembly/assembledGenome/Option_{outcome}_reads_basename.fasta"
    output:
        quastout="dir_to_each_assembly/assembledGenome/QUAST_OUT_{outcome}_reads_basename/report.txt"
    threads:
        threads_available
    shell:
        """
        echo "$(date): QUAST for reads_basename is starting!"
        """
        """
        quast.py -t {threads} -o dir_to_each_assembly/assembledGenome/QUAST_OUT_{wildcards.outcome}_reads_basename {input.option}
        """
        """
        echo "$(date): QUAST for reads_basename is finished!"
        """

rule nucmer:
    input:
        ref="path_to_ref_fasta",
        option="dir_to_each_assembly/assembledGenome/Option_{outcome}_reads_basename.fasta"
    output:
        delta="dir_to_each_assembly/prep/delta_{outcome}_reads_basename.delta",
        coords="dir_to_each_assembly/prep/coord_{outcome}_reads_basename.coords"
    threads:
        threads_available
    shell:
        """
        nucmer -t {threads} --delta={output.delta} {input.ref} {input.option}
        """
        """
        show-coords {output.delta} > {output.coords}
        """

rule samtools:
    input:
        option="dir_to_each_assembly/assembledGenome/Option_{outcome}_reads_basename.fasta",
        coords="dir_to_each_assembly/prep/coord_{outcome}_reads_basename.coords"
    output:
        region1="dir_to_each_assembly/prep/{outcome}_reads_basename_region1.fa",
        region2="dir_to_each_assembly/prep/{outcome}_reads_basename_region2.fa",
        full="dir_to_each_assembly/standardizedGenome/{outcome}_reads_basename_full.fa"
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
        full="dir_to_each_assembly/standardizedGenome/{outcome}_reads_basename_full.fa"
    output:
        outfiles="dir_to_each_assembly/annotatedGenome/{outcome}_reads_basename_full.gb"
    shell:
        """
        echo "$(date): PGA for reads_basename is starting!"
        """
        """
        perl path_to_PGA -r path_to_ref_gb -t dir_to_each_assembly/standardizedGenome -o dir_to_each_assembly/annotatedGenome
        """
        """
        echo "$(date): PGA for reads_basename is finished!"
        """
