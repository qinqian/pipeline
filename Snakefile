#input_sample = config['input'].split(',')
input_sample = "sample"

rule all:
   input:
       multiext("MN908947.3.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"),
       expand("mapped/{sample}.bam", sample=input_sample),
       expand("dedup/{sample}.bam", sample=input_sample),
       expand("dedup/{sample}.sorted.bam", sample=input_sample),
       expand("dedup/{sample}.sorted.bam.bai", sample=input_sample),
       expand("dedup/{sample}.metrics.txt", sample=input_sample),
       expand("samtools_stats/{sample}.txt", sample=input_sample),
       "qc/multiqc.html",
       expand("calls/{sample}.vcf", sample=input_sample)

rule bwa_index:
    input:
        "MN908947.3.fasta",
    output:
        idx=multiext("MN908947.3.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa_index/MN908947.3.log",
    params:
        algorithm="is",
    wrapper:
        "v1.1.0/bio/bwa/index"

rule samtools_fai_index:
    input:
        "MN908947.3.fasta"
    output:
        "MN908947.3.fasta.fai"
    wrapper:
        "v1.1.0/bio/samtools/faidx"

rule bwa_mem:
    input:
        reads=["{sample}.R1.paired.fq.gz", "{sample}.R2.paired.fq.gz"],
        idx=multiext("MN908947.3.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai")
    output:
        "mapped/{sample}.bam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sort="samtools",
        sort_order="queryname"
    threads: 8
    wrapper:
        "v1.1.0/bio/bwa/mem"

rule samtools_stats:
    input:
        "mapped/{sample}.bam"
    output:
        "samtools_stats/{sample}.txt"
    log:
        "logs/samtools_stats/{sample}.log"
    wrapper:
        "v1.1.0/bio/samtools/stats"

rule multiqc:
    input:
        expand("samtools_stats/{sample}.txt", sample=input_sample)
    output:
        "qc/multiqc.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "v1.1.0/bio/multiqc"

rule mark_duplicates:
    input:
        "mapped/{sample}.bam"
    output:
        bam="dedup/{sample}.bam",
        metrics="dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log",
    params:
        extra=" --REMOVE_DUPLICATES true --ASSUME_SORT_ORDER queryname ",
    resources:
        mem_mb=1024,
    wrapper:
        "v1.1.0/bio/picard/markduplicates"

rule samtools_sort:
    input:
        "dedup/{sample}.bam"
    output:
        "dedup/{sample}.sorted.bam"
    params:
        extra = "-m 8G",
        tmp_dir = "/tmp/"
    threads:
        8
    wrapper:
        "v1.1.0/bio/samtools/sort"

rule samtools_index:
    input:
        "dedup/{sample}.sorted.bam"
    output:
        "dedup/{sample}.sorted.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    wrapper:
        "v1.1.0/bio/samtools/index"

rule freebayes:
    input:
        ref="MN908947.3.fasta",
        idx=multiext("MN908947.3.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"),
        samples="dedup/{sample}.sorted.bam",
        indexes="dedup/{sample}.sorted.bam.bai"
    output:
        "calls/{sample}.vcf"
    log:
        "logs/freebayes/{sample}.log",
    params:
        extra="",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 2
    wrapper:
        "v1.1.0/bio/freebayes"

