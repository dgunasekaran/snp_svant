import os

rule sort_bam:
    input:
        os.path.join(config["outdir"], "preprocessed/mapped/{sample}_aligned.bam"),
    output:
        os.path.join(config["outdir"], "preprocessed/mapped/{sample}_aligned_sorted.bam"),
    log:
        os.path.join(config["logs"], "samtools_sort_{sample}.log"),
    params:
        extra="-m 4G",
    threads:
        config["threads"]
    wrapper:
        "v2.2.1/bio/samtools/sort"
