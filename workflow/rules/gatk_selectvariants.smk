import os

rule select_snps_round_1:
    input:
        vcf=os.path.join(config["outdir"], "preprocessed/3_hapcaller_r1/{sample}_raw_variants.vcf"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"], "preprocessed/3_hapcaller_r1/{sample}_raw_snps.vcf"),
    log:
        os.path.join(config["logs"], "gatk_select_snps_round1_{sample}.log"),
    params:
        extra="--select-type-to-include SNP",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/selectvariants"

rule select_indels_round_1:
    input:
        vcf = os.path.join(config["outdir"],"preprocessed/3_hapcaller_r1/{sample}_raw_variants.vcf"),
        ref = config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"],"preprocessed/3_hapcaller_r1/{sample}_raw_indels.vcf"),
    log:
        os.path.join(config["logs"],"gatk_select_indels_round1_{sample}.log"),
    params:
        extra="--select-type-to-include INDEL",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/selectvariants"

rule exclude_snps_round_1:
    input:
        vcf=os.path.join(config["outdir"], "preprocessed/3_hapcaller_r1/{sample}_filtered_snps.vcf"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"], "preprocessed/3_hapcaller_r1/{sample}_bqsr_snps.vcf"),
    log:
        os.path.join(config["logs"], "gatk_exclude_snps_round1_{sample}.log"),
    params:
        extra="--exclude-filtered",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/selectvariants"

rule exclude_indels_round_1:
    input:
        vcf=os.path.join(config["outdir"], "preprocessed/3_hapcaller_r1/{sample}_filtered_indels.vcf"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"], "preprocessed/3_hapcaller_r1/{sample}_bqsr_indels.vcf"),
    log:
        os.path.join(config["logs"], "gatk_exclude_indels_round1_{sample}.log"),
    params:
        extra="--exclude-filtered",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/selectvariants"

rule select_snps_round_2:
    input:
        vcf=os.path.join(config["outdir"],"preprocessed/hapcaller_r2/{sample}_raw_variants_recal.vcf"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"],"preprocessed/hapcaller_r2/{sample}_raw_snps_recal.vcf"),
    log:
        os.path.join(config["logs"],"gatk_select_snps_round2_{sample}.log"),
    params:
        extra="--select-type-to-include SNP",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/selectvariants"

rule select_indels_round_2:
    input:
        vcf=os.path.join(config["outdir"],"preprocessed/hapcaller_r2/{sample}_raw_variants_recal.vcf"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"],"preprocessed/hapcaller_r2/{sample}_raw_indels_recal.vcf"),
    log:
        os.path.join(config["logs"],"gatk_select_indels_round2_{sample}.log"),
    params:
        extra="--select-type-to-include INDEL",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/selectvariants"

rule exclude_snps_round_2:
    input:
        vcf=os.path.join(config["outdir"], "preprocessed/hapcaller_r2/{sample}_filtered_snps_final.vcf"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"], "preprocessed/final_variants/{sample}_filtered_snps_final.vcf"),
    log:
        os.path.join(config["logs"], "gatk_exclude_snps_round2_{sample}.log"),
    params:
        extra="--exclude-filtered",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/selectvariants"

rule exclude_indels_round_2:
    input:
        vcf=os.path.join(config["outdir"], "preprocessed/hapcaller_r2/{sample}_filtered_indels_final.vcf"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(config["outdir"], "preprocessed/final_variants/{sample}_filtered_indels_final.vcf"),
    log:
        os.path.join(config["logs"], "gatk_exclude_indels_round2_{sample}.log"),
    params:
        extra="--exclude-filtered",
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/gatk/selectvariants"
