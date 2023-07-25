import os

rule insert_size:
    input:
        os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_sorted_dedup_reads.bam"),
    output:
        txt=os.path.join(config["outdir"], "preprocessed/metrics/{sample}_insert_metrics.txt"),
        pdf=os.path.join(config["outdir"], "preprocessed/metrics/{sample}_insert_size_histgram.pdf"),
    log:
        os.path.join(config["logs"], "picard_insertsize_{sample}.log"),
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/picard/collectinsertsizemetrics"
