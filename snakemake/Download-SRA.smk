configfile: "/home/yanjun/prontfs/yanjunzan/Nutstore/home_bin/python/yaml/Download-SRA.yaml"

SAMPLES = ["ERR274529","ERR274543","ERR274526"]

WORKING_DIR = config['wd']
RESULT_DIR = config["rd"]
RESULT_DIR2 = config["rd"]

REF_dir = config["refdir"]
genome = config["ref"]


rule all:
    input:    
        #fq1 = expand(WORKING_DIR + "fastq/{sample}_1.fastq.gz",sample = SAMPLES),
        #fq2 = expand(WORKING_DIR + "fastq/{sample}_2.fastq.gz",sample = SAMPLES),
        bam = expand(RESULT_DIR2 + "mapped/{sample}_sort_index.bam",sample = SAMPLES)
        #cov = expand(RESULT_DIR2 + "cov/{sample}_cov.txt",sample = SAMPLES)
    message:
        "Job done! Removing temporary directory"

    
rule get_SRR_files:
    output:
        fw = temp(WORKING_DIR + "fastq/{sample}_1.fastq.gz"),
        rev= temp(WORKING_DIR + "fastq/{sample}_2.fastq.gz")
    params:
       SRA = "{sample}",
       DIR = config["wd"]+"fastq/",
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
        # 2> {log}
rule bwa_map:
    input:
        ref = REF_dir + genome + ".fa",
        fq1  = WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",
        fq2  = WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz"
    params:
        threads = config["threads"],
        ref = REF_dir + genome + ".fa"
        #bwamem2 = config["bwamem2"]
    output:
        #sort_index_bam = WORKING_DIR + "mapped/{sample}_sort_index.bam",
        sort_bam = RESULT_DIR2 + "mapped/{sample}_sort_index.bam"
    shell:
        "bwa  mem -t {params.threads} {params.ref} {input.fq1} {input.fq2} | \
        samtools sort -O BAM -@ {params.threads}  -o {output.sort_bam} "
        " && samtools index -@ {params.threads} {output.sort_bam}"
rule GetDep:
    input:
        sort_bam = RESULT_DIR2 + "mapped/{sample}_sort_index.bam"
    output:
        cov = RESULT_DIR2 + "cov/{sample}_cov.txt"
    shell:
        "samtools depth -a -q 30 -Q 30 {input} >> {output}" # filter 0,  | awk '$3 >0 {print }' | awk 'NR%10==0 {print $0 }' >> "

#rule qbam:
#    input:
#        sort_bam = RESULT_DIR2 + "mapped/{sample}_sort_index.bam"
#    output:
#        outfolder =
#        out
#    params:
#    shell:
#        "qualimap bamqc -bam ERR274526_sort_index.bam -c -outdir ./ "