import os

rule picard_sort:
    input:
        os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_picard_marked_duplicates.bam"),
    output:
        os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_sorted_dedup_reads.bam"),
    log:
        os.path.join(config["logs"], "picard_sort_{sample}.log"),
    params:
        sort_order="coordinate",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/picard/sortsam"
