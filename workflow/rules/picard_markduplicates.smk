import os

rule mark_duplicates:
    input:
        bams=os.path.join(config["outdir"], "preprocessed/mapped/{sample}_aligned_sorted.bam"),
    output:
        bam=os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_picard_marked_duplicates_wo_rg.bam"),
        metrics=os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_picard_marked_dup_metrics.txt"),
    log:
        os.path.join(config["logs"], "markdups_{sample}.log"),
    params:
        extra="--REMOVE_DUPLICATES false",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/picard/markduplicates"
