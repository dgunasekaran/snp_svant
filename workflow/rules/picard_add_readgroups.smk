import os

rule add_readgroups:
    input:
        os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_picard_marked_duplicates_wo_rg.bam"),
    output:
        os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_picard_marked_duplicates.bam"),
    log:
        os.path.join(config["logs"], "add_readgroups_{sample}.log"),
    params:
        extra="--RGID 1 --RGLB library1 --RGPL illumina --RGPU {sample} --RGSM {sample}",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/picard/addorreplacereadgroups"
