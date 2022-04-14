rule all:
    input:
        "WORKDIR/annotatedGenome/prefix_check_stop_codons.txt"
        
rule decompress:
    input:
        compressed_forward="read1", 
        compressed_reverse="read2"
    output:
        decompressed_forward="WORKDIR/reads/forward_read.fastq",
        decompressed_reverse="WORKDIR/reads/reverse_read.fastq"
    threads:
        threads_available
    shell:
        """
        echo "$(date): Starts decompressing prefix!"
        pigz -d -p {threads} -k -c {input.compressed_forward} > {output.decompressed_forward}
        pigz -d -p {threads} -k -c {input.compressed_reverse} > {output.decompressed_reverse}
        echo "$(date): Finished decompressing prefix!"
        """

rule trimmomatic:
    input:
        forward_read="WORKDIR/reads/forward_read.fastq",
        reverse_read="WORKDIR/reads/reverse_read.fastq"
    output:
        forwardPaired="WORKDIR/trimmedReads/forward_readP.fastq",
        reversePaired="WORKDIR/trimmedReads/reverse_readP.fastq", 
        forwardUnPaired="WORKDIR/trimmedReads/forward_readU.fastq", 
        reverseUnPaired="WORKDIR/trimmedReads/reverse_readU.fastq"
    threads:
        threads_available
    log: "WORKDIR/logs/trimmomatic.log"
    shell:
        """
        echo "$(date): Trimmomatic for prefix is starting!"
        java -jar path_to_trimmomatic PE -threads {threads} {input.forward_read} {input.reverse_read} {output.forwardPaired} {output.forwardUnPaired} {output.reversePaired} {output.reverseUnPaired} ILLUMINACLIP:ADAPTERS:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60 &> {log}
        echo "$(date): Trimmomatic for prefix is finished!"
        """

rule novoplasty:
    input:
        config_novo="WORKDIR/prefix_config.txt",
        fowardPaired="WORKDIR/trimmedReads/forward_readP.fastq",
        reversePaired="WORKDIR/trimmedReads/reverse_readP.fastq"
    output:
        directory("WORKDIR/assembledGenome")
    log: "WORKDIR/logs/novoplasty.log"
    shell:
        """
        echo "$(date): NOVOPlasty for prefix is starting!"
        mkdir {output}
        perl path_to_novoplasty -c {input.config_novo} &> {log}
        echo "$(date): NOVOPlasty for prefix is finished!"
        """

rule standardize:
    input:
        "WORKDIR/assembledGenome"
    output:
        directory("WORKDIR/standardizedGenome"),
        standardized="WORKDIR/standardizedGenome/prefix.plastome_assembly.fa",
    log: "WORKDIR/logs/standardize.log"
    shell:
        """
        path_to_repo/standardize_cpDNA.sh -d {input} -o {output.standardized} -p prefix &> {log}
        """

rule PGA:
    input:
        "WORKDIR/standardizedGenome"
    output:
        directory("WORKDIR/annotatedGenome"),
        gb="WORKDIR/annotatedGenome/test.plastome_assembly.gb"
    log: "WORKDIR/logs/PGA.log"
    shell:
        """
        echo "$(date): PGA for prefix is starting!"
        perl path_to_PGA -r path_to_ref_gb -t {input} -o {output} &> {log}
        echo "$(date): PGA for prefix is finished!"
        """
      
rule ISC:
    input:
        fa="WORKDIR/standardizedGenome/prefix.plastome_assembly.fa",
        gb="WORKDIR/annotatedGenome/prefix.plastome_assembly.gb"
    output:
        out="WORKDIR/annotatedGenome/prefix_check_stop_codons.txt"
    shell:
        """
        echo "$(date): Checking for internal stop codons"
        WORKDIR/prefix_isc.sh {input.gb} {input.fa} > {output.out}
        echo "$(date): DONE: check the output in WORKDIR/annotatedGenome/"
        """

