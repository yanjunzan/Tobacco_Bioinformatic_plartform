configfile: "/home/yanjun/prontfs/yanjunzan/Nutstore/home_bin/python/yaml/RNAseq.yaml"

import pandas as pd

sra = pd.read_csv("/home/yanjun/prontfs/project/tobacco/results/SraRunTable2.txt",sep="\t")
SAMPLES = sra["Run"].tolist()#[801:1200]

WORKING_DIR = config['working_dir']
RESULT_DIR = config["working_dir"]
REF_dir = config["refdir"]
genome = config["ref"]


rule all:
    input:    
        #fq1 = expand(WORKING_DIR + "fastq/{sample}_1.fastq.gz",sample = SAMPLES),
        #fq2 = expand(WORKING_DIR + "fastq/{sample}_2.fastq.gz",sample = SAMPLES),
        #bam = expand(RESULT_DIR + "mapped/{sample}.bam",sample = SAMPLES),
        #txt = expand(WORKING_DIR + "Expression/{sample}.txt",sample = SAMPLES),
        tsv = expand(WORKING_DIR + "tsv/{sample}.tsv",sample = SAMPLES)
    message:
        "Job done! Removing temporary directory"

    
rule get_SRR_files:
    output:
        fw = temp(WORKING_DIR + "fastq/{sample}_1.fastq.gz"),
        rev= temp(WORKING_DIR + "fastq/{sample}_2.fastq.gz")
    params:
       SRA = "{sample}",
       DIR = WORKING_DIR+"fastq/",
       tmpd = "/home/yanjun/prontfs/"
    message:
        "using fastq-dump to download SRA data files to {output.fw}"
    shell:
        "parallel-fastq-dump -T {params.tmpd} -t 10  -O {params.DIR} --split-files --gzip -s {params.SRA}"

rule fastp:
    input:
        fw = WORKING_DIR + "fastq/{sample}_1.fastq.gz",
        rev= WORKING_DIR + "fastq/{sample}_2.fastq.gz"
    output:
        fq1  = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",
        html = RESULT_DIR + "fastp/{sample}.html"
    message:"trimming {wildcards.sample} reads to {output.fq1}"
    threads: config["threads"]
    log:
        WORKING_DIR + "logs/fastp/{sample}.log.txt"
    params:
        qualified_quality_phred = config["fastp"]["qualified_quality_phred"]
    shell:
        "fastp --thread {threads}  --html {output.html} \
        --qualified_quality_phred {params.qualified_quality_phred} \
        --in1 {input.fw} --in2 {input.rev} --out1 {output.fq1} --out2 {output.fq2} \
        > {log} 2>&1"

rule hisat_mapping:
    input:
        fq1  = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    output:
        bams  = WORKING_DIR + "mapped/{sample}.bam",
        met   = RESULT_DIR + "logs/{sample}_met.txt"
    params:
        indexName = REF_dir + genome,
        threads = config["threads"]
    message:
        "mapping reads to genome to bam files."
    shell:
        "hisat2 -p {threads} --met-file {output.met} -x {params.indexName} \
            -1 {input.fq1} -2 {input.fq2} | samtools sort -O BAM -@ {threads}  -o {output.bams} "
            "&& samtools index -@ {params.threads} {output.bams}"

rule strigtie:
    input:
        GFF =  config["gff"],
        bams  = WORKING_DIR + "mapped/{sample}.bam"
    output:
        gif =  WORKING_DIR + "gif/{sample}.gif",
        tsv =  WORKING_DIR + "tsv/{sample}.tsv"
    shell:
        "stringtie -p 5 -e -B -G {input.GFF} \
        -o {output.gif} \
        -A {output.tsv} {input.bams}"

#rule featurecount:
#    input:
#        GFF =  config["gff"],
#        bams  = WORKING_DIR + "mapped/{sample}.bam"
#    output:
#        txt= WORKING_DIR + "Expression/{sample}.txt"
#
#    shell:
#        "featureCounts -T 10 -p -t exon -g transcript_id -a {input.GFF} \
#        -o  {output.txt} {input.bams} "