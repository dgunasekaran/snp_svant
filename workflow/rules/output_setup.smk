import os

rule fq_dump:
    output:
        base_dir=os.path.join(config["outdir"],"preprocessed/"),
        align_dir=os.path.join(config["outdir"],"preprocessed/mapped"),
        dup_dir=os.path.join(config["outdir"],"preprocessed/markduplicates"),
        metric_dir=os.path.join(config["outdir"],"preprocessed/metrics"),
        hap_dir_1=os.path.join(config["outdir"],"preprocessed/3_hapcaller_r1"),
        bqsr_dir_1=os.path.join(config["outdir"],"preprocessed/bqsr_r1"),
        hap_dir_2=os.path.join(config["outdir"],"preprocessed/hapcaller_r2"),
        bqsr_dir_2=os.path.join(config["outdir"],"preprocessed/bqsr_r2"),
        variant_dir=os.path.join(config["outdir"],"preprocessed/final_variants"),
        gridss_dir=os.path.join(config["outdir"],"preprocessed/gridss"),
        vep_dir=os.path.join(config["outdir"],"preprocessed/vep_genes"),
    params:
        out_dir=config["outdir"],
    shell:
        '''
        mkdir -p {params.out_dir}
        mkdir -p {config[logs]}
        mkdir -p {output.base_dir}
        mkdir -p {output.align_dir}
        mkdir -p {output.dup_dir}
        mkdir -p {output.metric_dir}
        mkdir -p {output.hap_dir_1}
        mkdir -p {output.hap_dir_2}
        mkdir -p {output.bqsr_dir_1}
        mkdir -p {output.bqsr_dir_2}
        mkdir -p {output.variant_dir}
        mkdir -p {output.vep_dir}
        mkdir -p {output.gridss_dir}
        '''

