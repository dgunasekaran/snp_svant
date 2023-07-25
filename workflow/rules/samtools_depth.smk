import os

rule depth_summary:
    input:
        bams=[os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_sorted_dedup_reads.bam")],
    output:
        os.path.join(config["outdir"], "preprocessed/metrics/{sample}_depth_out.txt"),
    log:
        os.path.join(config["logs"],"samtools_depth_{sample}.log"),
    wrapper:
        "v2.2.1/bio/samtools/depth"
