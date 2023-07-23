import os

if config["import"]:
    rule fq_dump:
        output:
            os.path.join(config["outdir"], "raw/{sample}_1.fastq"),
            os.path.join(config["outdir"], "raw/{sample}_2.fastq")
        log:
            os.path.join(config["logs"], "data_dump_{sample}.log")
        params:
            extra="--skip-technical"
        threads:
            config["threads"]
        wrapper:
            "v2.2.1/bio/sra-tools/fasterq-dump"
