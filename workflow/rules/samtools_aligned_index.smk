import os

rule samtools_aligned_index:
    input:
        os.path.join(config["outdir"], "preprocessed/mapped/{sample}_aligned_sorted.bam"),
    output:
        os.path.join(config["outdir"], "preprocessed/mapped/{sample}_aligned_sorted.bam.bai"),
    log:
        os.path.join(config["logs"], "samtools_aligned_index_{sample}.log"),
    threads: config["threads"]
    wrapper:
        "v2.2.1/bio/samtools/index"
