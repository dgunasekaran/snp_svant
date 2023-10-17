import os


rule vep_custom:
    input:
        vcf=os.path.join(config["outdir"], "preprocessed/final_variants/{sample}_filtered_snps_final.vcf"),
    output:
        vcf_gz=os.path.join(config["outdir"], "preprocessed/final_variants/{sample}_filtered_snps_final.vcf.gz"),
        stats_file=os.path.join(config["outdir"],"preprocessed/vep_genes/{sample}_snps.html"),
        vep_txt=os.path.join(config["outdir"],"preprocessed/vep_genes/{sample}_snps_vep.txt"),
    params:
        gff_gz=config["reference"]["gff"]+".gz",
        fasta_gz=config["reference"]["genome"]+".gz",
    log:
        os.path.join(config["logs"],"vep_{sample}.log")
    threads:
        config["threads"]
    shell:
        '''
        grep -v "#" {config[reference][gff]} | sort -k1,1 -k4,4n | bgzip -c > {params.gff_gz}
        tabix -p gff {params.gff_gz}
        bgzip -c {config[reference][genome]} > {params.fasta_gz}
        bgzip -c {input.vcf} > {output.vcf_gz}
        echo {wildcards.sample}
        ./external_tools/ensembl-vep/vep -i {output.vcf_gz} --fasta {params.fasta_gz} --gff {params.gff_gz} --stats_file {output.stats_file} -o {output.vep_txt}
        '''

