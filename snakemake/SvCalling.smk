infolder = "/home/yanjun/proext4/tobacco/results/SV/in/"
outfolder = "/home/yanjun/proext4/tobacco/results/SV/out/"
ref = "ZY300.fa.gz"
reffai = "ZY300.fa.gz.fai"

import os,re
a =  os.listdir(infolder)

Samples = [re.sub("(.*)_sort_index_rg_nodup\\.bam","\\1",f) for f in a if re.match(".*nodup\\.bam$",f) ]

rule all:
    input:
         #xpand(WORKING_DIR + "fastq/{sample}_2.fastq.gz",sample = SAMPLES),
        output = expand( outfolder +"{sample}"+ "/combined.genotyped.vcf",sample = Samples)
    message:
        "Job done! Removing temporary directory"

rule SvCall:
    input:
        bam =infolder + "{sample}"+ "_sort_index_rg_nodup.bam",
        
    output:
        outvcf = outfolder + "{sample}" + "/combined.genotyped.vcf"
    params:
        inF = infolder,
        bam = "{sample}"+ "_sort_index_rg_nodup.bam",
        bamfai = "{sample}" + "_sort_index_rg_nodup.bai",
        ref = "ZY300.fa.gz",
        reffai = "ZY300.fa.gz.fai",
        infolder1=outfolder + "{sample}"
    message:
        "processing {input} for sv calling "
    shell:
        " mkdir -p {params.infolder1} && \
        docker run -v {params.inF}:/home/dnanexus/in -v {params.infolder1}:/home/dnanexus/out \
        dnanexus/parliament2 \
        --bam {params.bam} \
        --bai {params.bamfai} \
        --fai {params.reffai} \
        -r {params.ref} \
       --delly_duplication --delly_inversion --delly_insertion --delly_deletion \
       --cnvnator --manta --genotype "
