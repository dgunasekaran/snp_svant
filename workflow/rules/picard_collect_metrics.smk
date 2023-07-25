import os

rule alignment_summary:
    input:
        ref=config["reference"]["genome"],
        bam=os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_sorted_dedup_reads.bam"),
    output:
        os.path.join(config["outdir"], "preprocessed/metrics/{sample}_alignment_metrics.txt"),
    log:
        os.path.join(config["logs"], "picard_summary_{sample}.log"),
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/picard/collectalignmentsummarymetrics"
