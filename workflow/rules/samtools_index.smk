import os

rule samtools_index:
    input:
        os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_sorted_dedup_reads.bam"),
    output:
        os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_sorted_dedup_reads.bam.bai"),
    log:
        os.path.join(config["logs"], "samtools_index_{sample}.log"),
    threads: config["threads"]
    wrapper:
        "v2.2.1/bio/samtools/index"
