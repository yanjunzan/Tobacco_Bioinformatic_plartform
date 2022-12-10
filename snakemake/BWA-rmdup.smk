configfile: "./yaml/bwamap.yaml"


WORKING_DIR = config["working_dir"]
RESULT_DIR = config["result_dir"]
RESULT_DIR2 = config["result_dir2"]

REF_dir = config["refdir"]
genome = config["ref"]

SAMPLES = [re.sub("(.*)_R.*","\\1",f) for f in os.listdir(WORKING_DIR+"trimmed") if re.search("1r_R",f)]



rule all:
    input:
        #vcf = WORKING_DIR + "vcf/all.vcf",
        #reffai = REF_dir + genome + ".fasta.fai",
        #bamlist= WORKING_DIR + "mapped/BAMList.txt",
        outbam = expand(RESULT_DIR2 + "mapped/{sample}_sort_index_rg_nodup.bam",sample=SAMPLES),
        outbai = expand(RESULT_DIR2 + "mapped/{sample}_sort_index_rg_nodup.bam.bai",sample=SAMPLES)

        #samstats = expand(WORKING_DIR + "bamstats/{sample}.bamstats.stat",sample=SAMPLES),
        #sort_index_rgbam = expand(WORKING_DIR + "mapped/{sample}_sort_index_rg.bam",sample=SAMPLES),
        #genomeF = REF_dir + genome + ".bwt",
        #qcfile = expand(RESULT_DIR + "fastp/{sample}.html",sample=SAMPLES),
        #fq1 = expand(WORKING_DIR + "trimmed/{sample}_R1_trimmed.fq.gz",sample = SAMPLES),
        #fq2 = expand(WORKING_DIR + "trimmed/{sample}_R2_trimmed.fq.gz",sample = SAMPLES),
        #bam = expand(WORKING_DIR + "mapped/{sample}_sort_index.bam", sample = SAMPLES)
    message:
        "Job done! Removing temporary directory"

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
        sort_bam = temp(RESULT_DIR2 + "mapped/{sample}_sort_index.bam"),
        tmpbam  = temp(WORKING_DIR + "trimmed/{sample}_trimmed.fq.gz")
    shell:
        "cat {input.fq1} {input.fq2} > {output.tmpbam}"
        " && bwa  mem -t {params.threads} {params.ref}  {output.tmpbam} | \
        samtools sort -O BAM -@ {params.threads}  -o {output.sort_bam} "
        " && samtools index -@ {params.threads} {output.sort_bam}"
rule addReadGroup:
    input:
        picard = config["picard"],
        sort_index_bam = RESULT_DIR2 + "mapped/{sample}_sort_index.bam"
    output:
        sort_index_rgbam = temp(RESULT_DIR2 + "mapped/{sample}_sort_index_rg.bam")
    params:
        sample = "{sample}"
    shell:
        "java -jar {input.picard} AddOrReplaceReadGroups \
        I={input.sort_index_bam}  O={output.sort_index_rgbam} RGID={params.sample} RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM={params.sample}"
#        
#rule samtools_stat:
#    input:
#        bam = RESULT_DIR2 + "mapped/{sample}_sort_index_rg.bam"
#    output:
#        samstats = RESULT_DIR2 + "bamstats/{sample}.bamstats.stat"
#    params:
#        threads= config["threads"]
#    shell:
#        "samtools stats -@ {params.threads} {input.bam} > {output.samstats}"
#        
rule picard_remove_duplicates:
    input:
        picard = config["picard"],
        bam = RESULT_DIR2 + "mapped/{sample}_sort_index_rg.bam"
    output:
        outbam = protected(RESULT_DIR2 + "mapped/{sample}_sort_index_rg_nodup.bam"),
        metrics = RESULT_DIR2 + "dupstats/{sample}.picard.marked_dup_metrics.txt",
        outbai = protected(RESULT_DIR2 + "mapped/{sample}_sort_index_rg_nodup.bai")
    shell:
        "java -Xmx6G -jar {input.picard}  MarkDuplicates -I {input.bam} -O {output.outbam} -M {output.metrics} --REMOVE_DUPLICATES true --CREATE_INDEX true"
        
        
rule fix_idx:
    input:
        outbai = RESULT_DIR2 + "mapped/{sample}_sort_index_rg_nodup.bai",

    output:
        outbai = protected(RESULT_DIR2 + "mapped/{sample}_sort_index_rg_nodup.bam.bai"),
    shell:
        "mv {input.outbai} {output.outbai}"
        