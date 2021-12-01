wildcard_constraints:
    outcome="\d+"

OUTCOME= ["1","2"]

rule all:
    input:
        expand("WORKDIR/annotatedGenome/Option_{outcome}_check_stop_codons.txt", outcome=OUTCOME)
        
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
        """
        """
        pigz -d -p {threads} -k -c {input.compressed_forward} > {output.decompressed_forward}
        """
        """
        pigz -d -p {threads} -k -c {input.compressed_reverse} > {output.decompressed_reverse}
        """
        """
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
        """
        """
        java -jar path_to_trimmomatic PE -threads {threads} {input.forward_read} {input.reverse_read} {output.forwardPaired} {output.forwardUnPaired} {output.reversePaired} {output.reverseUnPaired} ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:60 &> {log}
        """
        """
        echo "$(date): Trimmomatic for prefix is finished!"
        """

rule novoplasty:
    input:
        config_novo="WORKDIR/prefix_config.txt",
        fowardPaired="WORKDIR/trimmedReads/forward_readP.fastq",
        reversePaired="WORKDIR/trimmedReads/reverse_readP.fastq"
    output:
        option1="WORKDIR/assembledGenome/Option_1_prefix.fasta",
        option2="WORKDIR/assembledGenome/Option_2_prefix.fasta"
    log: "WORKDIR/logs/novoplasty.log"
    shell:
        """
        echo "$(date): NOVOPlasty for prefix is starting!"
        """
        """
        perl path_to_novoplasty -c {input.config_novo} &> {log}
        """
        """
        echo "$(date): NOVOPlasty for prefix is finished!"
        """

rule standardize:
    input:
        option="WORKDIR/assembledGenome/Option_{outcome}_prefix.fasta"
    output:
        standardized="WORKDIR/standardizedGenome/Option_{outcome}_prefix.standardized.fa",
    shell:
        """
        path_to_repo/standardize_cpDNA.sh -i {input.option} -o ${output.standardized}
        """

rule PGA:
    input:
        full="WORKDIR/standardizedGenome/Option_{outcome}_prefix.standardized.fa"
    output:
        outfiles="WORKDIR/annotatedGenome/Option_{outcome}_prefix_full.gb"
    log: "WORKDIR/logs/PGA_{outcome}.log"
    shell:
        """
        echo "$(date): PGA for prefix is starting!"
        """
        """
        perl path_to_PGA -r path_to_ref_gb -t WORKDIR/standardizedGenome/ -o WORKDIR/annotatedGenome &> {log}
        """
        """
        echo "$(date): PGA for prefix is finished!"
        """
      
rule ISC:
    input:
        fa="WORKDIR/standardizedGenome/Option_{outcome}_prefix.standardized.fa",
        gb="WORKDIR/annotatedGenome/Option_{outcome}_prefix_full.gb"
    output:
        out="WORKDIR/annotatedGenome/Option_{outcome}_check_stop_codons.txt"
    shell:
        """
        path_to_repo/check_internal_stops.sh {input.gb} {input.fa} > {output.out}
        """
