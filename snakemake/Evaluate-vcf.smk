vcf = "GBS_5558_bi-SNP-dp2-miss.5-maf.03-fixdip-IMP-nowild.recode.vcf.gz"
suffix = "impute_all-nowild"
WORKING_DIR= "/home/yanjun/proext4/tobacco/results/"

rule all:
    input:
        imiss = WORKING_DIR + "vcf/" + suffix + ".imiss",
        lmiss = WORKING_DIR + "vcf/" + suffix + ".lmiss",
        freq = WORKING_DIR + "vcf/" + suffix + ".frq",
        snpden = WORKING_DIR + "vcf/" + suffix + ".snpden",
        het = WORKING_DIR + "vcf/" + suffix + ".het"
        #bin = WORKING_DIR + "plink.grm.bin"

    message:
        "Job done! Removing temporary directory"

rule CallRate:
    input:
        vcffile = WORKING_DIR + vcf
    output:
        imiss = WORKING_DIR + "vcf/" + suffix + ".imiss",
    params:
        suf  = WORKING_DIR + "vcf/" + suffix
    message:
        "processing {input} for individual missing rate!"
    shell:
        "vcftools --gzvcf {input} --missing-indv --out {params.suf}"
        
rule CallRate2:
    input:
        vcffile = WORKING_DIR + vcf
    output:
        lmiss = WORKING_DIR + "vcf/" + suffix + ".lmiss",
    params:
        suf  = WORKING_DIR + "vcf/" + suffix
    message:
        "processing {input} for locus missing rate!"
    shell:
        "vcftools --gzvcf {input} --missing-site --out {params.suf}"
       
rule freq:
    input:
        vcffile = WORKING_DIR + vcf
    output:
        freq = WORKING_DIR + "vcf/" + suffix + ".frq",
    params:
        suf  = WORKING_DIR + "vcf/" + suffix
    message:
        "processing {input} for site freq"
    shell:
        "vcftools --gzvcf {input} --freq --out {params.suf}"

rule SNPden:
    input:
        vcffile = WORKING_DIR + vcf
    output:
        snpden = WORKING_DIR + "vcf/" + suffix + ".snpden"
    message:
        "processing {input} for snp density!"
    params:
        suf  = WORKING_DIR + "vcf/" + suffix

    shell:
        "vcftools --gzvcf {input} --SNPdensity 5000000 --out {params.suf}"

rule Het:
    input:
        vcffile = WORKING_DIR + vcf
    output:
        het = WORKING_DIR + "vcf/" + suffix + ".het",
    params:
        suf  = WORKING_DIR + "vcf/" + suffix
    message:
        "processing {input} for het!"
    shell:
        "vcftools --gzvcf {input} --het --out {params.suf}"
        
rule QC:
    input:
        vcffile = WORKING_DIR + vcf
    output:
        het = WORKING_DIR + "vcf/" + suffix + ".recode.vcf"
    params:
        suf  = WORKING_DIR + "vcf/" + suffix
    message:
        "QC {input} "
    shell:
        "vcftools --gzvcf {input} --remove-indels --min-alleles 2 --max-alleles 2 \
        --max-missing 0.5 --min-meanDP 2 --max-maf 0.05 --recode --out {params.suf}"

rule GRM:
    input:
        vcffile = WORKING_DIR + vcf
    output:
        bin = WORKING_DIR + "plink.grm.bin"
    message:
        "processing {input} for grm!"
    shell:
        "plink2 --vcf {input} --make-grm-bin  --allow-extra-chr "
