import os


rule download_vep_plugins:
    output:
        directory("resources/vep/plugins")
    params:
        release=100
    wrapper:
        "v2.6.0/bio/vep/plugins"


rule zip_vcf:
    input:
        vcf=os.path.join(config["outdir"],"preprocessed/final_variants/{sample}_filtered_snps_final.vcf"),
    output:
        vcf_gz=os.path.join(config["outdir"],"preprocessed/final_variants/{sample}_filtered_snps_final.vcf.gz"),
    shell:
        '''
        bgzip -c {input.vcf} > {output.vcf_gz}
        '''


rule vep_wrapper:
    input:
        calls=os.path.join(config["outdir"], "preprocessed/final_variants/{sample}_filtered_snps_final.vcf.gz"),
        fasta=config["reference"]["genome"],
        fai=config["reference"]["genome"]+".fai", # fasta index
        gff=config["reference"]["gff"],
        plugins="resources/vep/plugins",
    output:
        calls=os.path.join(config["outdir"],"preprocessed/vep_genes/{sample}_snps_vep.vcf"),
        stats=os.path.join(config["outdir"],"preprocessed/vep_genes/{sample}_snps.html"),
    params:
        plugins = ["LoFtool"],
        extra="--everything",  # optional: extra arguments
    log:
        os.path.join(config["logs"],"vep_{sample}.log")
    threads:
        4
    wrapper:
        "v2.6.0/bio/vep/annotate"