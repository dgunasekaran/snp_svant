import os


rule vep_wrapper:
    input:
        calls=os.path.join(config["outdir"], "preprocessed/final_variants/{sample}_filtered_snps_final.vcf"),
        fasta=config["reference"]["genome"],
        fai=config["reference"]["genome"]+".fai", # fasta index
        gff=config["reference"]["gff"],
    output:
        calls=os.path.join(config["outdir"],"preprocessed/vep_genes/{sample}_snps_vep.vcf"),
        stats=os.path.join(config["outdir"],"preprocessed/vep_genes/{sample}_snps.html"),
    params:
        extra="--everything",  # optional: extra arguments
    log:
        os.path.join(config["logs"],"vep_{sample}.log")
    threads:
        config["threads"]
    wrapper:
        "v2.6.0/bio/vep/annotate"