import os

rule bqsr_round_1:
    input:
        bam=os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_sorted_dedup_reads.bam"),
        ref=config["reference"]["genome"],
        dict=os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.dict'),
        known=[os.path.join(config["outdir"],"preprocessed/3_hapcaller_r1/{sample}_bqsr_snps.vcf"),
               os.path.join(config["outdir"],"preprocessed/3_hapcaller_r1/{sample}_bqsr_indels.vcf")],
    output:
        recal_table=os.path.join(config["outdir"],"preprocessed/bqsr_r1/{sample}_recal_data.table"),
    log:
        os.path.join(config["logs"],"gatk_bqsr_round1_{sample}.log"),
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/baserecalibrator"

rule apply_bqsr_round_1:
    input:
        bam=os.path.join(config["outdir"], "preprocessed/markduplicates/{sample}_sorted_dedup_reads.bam"),
        ref=config["reference"]["genome"],
        dict=os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.dict'),
        recal_table=os.path.join(config["outdir"],"preprocessed/bqsr_r1/{sample}_recal_data.table"),
    output:
        bam=os.path.join(config["outdir"],"preprocessed/bqsr_r1/{sample}_recal_reads.bam"),
    log:
        os.path.join(config["logs"],"gatk_apply_bqsr_round1_{sample}.log"),
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/applybqsr"

rule bqsr_round_2:
    input:
        bam=os.path.join(config["outdir"],"preprocessed/bqsr_r1/{sample}_recal_reads.bam"),
        ref=config["reference"]["genome"],
        dict=os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.dict'), \
        known = [os.path.join(config["outdir"],"preprocessed/3_hapcaller_r1/{sample}_bqsr_snps.vcf"),
                os.path.join(config["outdir"],"preprocessed/3_hapcaller_r1/{sample}_bqsr_indels.vcf")],
    output:
        recal_table=os.path.join(config["outdir"],"preprocessed/bqsr_r2/{sample}_post_recal_data.table"),
    log:
        os.path.join(config["logs"],"gatk_bqsr_round2_{sample}.log"),
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/baserecalibrator"
