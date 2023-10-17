import os


rule vep_custom:
    input:
        vcf=os.path.join(config["outdir"], "preprocessed/final_variants/{sample}_filtered_snps_final.vcf"),
    output:
        vcf_gz=os.path.join(config["outdir"], "preprocessed/final_variants/{sample}_filtered_snps_final.vcf.gz"),
        gff_gz=os.path.join(config["reference"]["gff"], ".gz"),
        fasta_gz=os.path.join(config["reference"]["genome"], ".gz"),
        stats_file=os.path.join(config["outdir"],"preprocessed/vep_genes/{sample}_snps.html"),
        vep_txt=os.path.join(config["outdir"],"preprocessed/vep_genes/{sample}_snps_vep.txt"),
    log:
        os.path.join(config["logs"],"vep_{sample}.log")
    threads:
        config["threads"]
    shell:
        '''
        grep -v "#" {config["reference"]["gff"]} | sort -k1,1 -k4,4n | bgzip -c > {output.gff_gz}
        tabix -p gff {output.gff_gz}
        bgzip -c {config["reference"]["genome"]} > {output.fasta_gz}
        bgzip -c {input.vcf} > {output.vcf_gz}
        echo {wildcards.sample}
        ./external_tools/ensembl-vep/vep -i {output.vcf_gz} --fasta {output.fasta_gz} --gff {output.gff_gz} --stats_file {output.stats_file} -o {output.vep_txt}
        '''

