import os

rule filter_snps_round_1:
    input:
        vcf=os.path.join(config["outdir"],"preprocessed/3_hapcaller_r1/{sample}_raw_snps.vcf"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"],"preprocessed/3_hapcaller_r1/{sample}_filtered_snps.vcf"),
    log:
        os.path.join(config["logs"],"gatk_filter_snps_round1_{sample}.log"),
    params:
        filters=config["snp_filters"],
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/variantfiltration"

rule filter_indels_round_1:
    input:
        vcf=os.path.join(config["outdir"],"preprocessed/3_hapcaller_r1/{sample}_raw_indels.vcf"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"],"preprocessed/3_hapcaller_r1/{sample}_filtered_indels.vcf"),
    log:
        os.path.join(config["logs"],"gatk_filter_indels_round1_{sample}.log"),
    params:
        filters=config["indel_filters"],
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/variantfiltration"

rule filter_snps_round_2:
    input:
        vcf=os.path.join(config["outdir"],"preprocessed/hapcaller_r2/{sample}_raw_snps_recal.vcf"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"],"preprocessed/hapcaller_r2/{sample}_filtered_snps_final.vcf"),
    log:
        os.path.join(config["logs"],"gatk_filter_snps_round2_{sample}.log"),
    params:
        filters=config["snp_filters"],
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/variantfiltration"

rule filter_indels_round_2:
    input:
        vcf=os.path.join(config["outdir"],"preprocessed/hapcaller_r2/{sample}_raw_indels_recal.vcf"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"],"preprocessed/hapcaller_r2/{sample}_filtered_indels_final.vcf"),
    log:
        os.path.join(config["logs"],"gatk_filter_indels_round2_{sample}.log"),
    params:
        filters=config["indel_filters"],
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/variantfiltration"
