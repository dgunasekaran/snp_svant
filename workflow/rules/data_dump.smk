import os

if config["import"]:
    rule fq_dump:
        output:
            os.path.join(config["outdir"],"raw/{sample}_1.fastq"),
            os.path.join(config["outdir"],"raw/{sample}_2.fastq")
        params:
            out_dir=config["outdir"],
            data_dir=config["outdir"] + "raw/",
        log:
            os.path.join(config["logs"],"data_dump_{sample}.log")
        threads:
            config["threads"]
        conda:
            "snp_svaunt"
        shell:
            '''
            mkdir -p {params.out_dir}
            mkdir -p {params.data_dir}
            mkdir -p {config[logs]}
            echo {wildcards.sample}
            ./external_tools/sratoolkit.3.0.6-mac64/bin/fasterq-dump {wildcards.sample} -e {threads} 2>&1 {log}
            mv {wildcards.sample}*.fastq {params.data_dir}
            '''

