import os


#reference_index_endings= (".amb",".ann",".bwt",".pac",".sa",".gridsscache",".img")
reference_index_endings= (".amb", ".ann", ".bwt", ".dict", ".fai", ".pac", ".sa")

rule gridss_custom:
    input:
        align_bam = os.path.join(config["outdir"],"preprocessed/mapped/{sample}_aligned_sorted.bam"),
        align_idx = os.path.join(config["outdir"],"preprocessed/mapped/{sample}_aligned_sorted.bam.bai"),
        reference = config["reference"]["genome_gridss"],
        dictionary = config["reference"]["genome_gridss"] + ".dict",
        refindex = multiext(config["reference"]["genome_gridss"],".amb",".ann",".bwt",".pac",".sa"),
    output:
        bam = os.path.join(config["outdir"],"preprocessed/gridss/{sample}_aligned.bam"),
        vcf=os.path.join(config["outdir"], "preprocessed/gridss/{sample}.vcf"),
        vep=os.path.join(config["outdir"], "preprocessed/gridss/{sample}"),
    params:
        out_dir=os.path.join(config["outdir"], "preprocessed/gridss/")
    log:
        os.path.join(config["logs"],"gridss_{sample}.log")
    threads:
        config["threads"]
    shell:
        '''
        mkdir -p {params.out_dir}
        mkdir -p {config[logs]}
        echo {wildcards.sample}
        ./external_tools/gridss_v_2_12_0/gridss -j external_tools/gridss_v_2_12_0/gridss-2.12.0-gridss-jar-with-dependencies.jar -r {input.reference} \
        -o {output.vcf} -a {output.bam} -w {params.out_dir} {input.align_bam}
        Rscript wworkflow/scripts/sv_annotation.R --vcf {output.vcf} --output_prefix {output.vep}
        '''

