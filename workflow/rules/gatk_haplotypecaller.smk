import os

rule variant_calling_round_1:
    input:
        bam=os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_sorted_dedup_reads.bam"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"], "preprocessed/hapcaller_r1/{sample}_raw_variants.vcf"),
    log:
        os.path.join(config["logs"], "gatk_hc_round1_{sample}.log"),
    threads: config["threads"]
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/haplotypecaller"
