import os
import sys

# Map system → binary path
SRA_BINARIES = {
    "macOS": "./external_tools/sratoolkit.3.0.6-mac64/bin/fasterq-dump",
    "linux": "./external_tools/sratoolkit.3.2.1-ubuntu64/bin/fasterq-dump",
}

if config["import"]:

    try:
        FASTERQ = SRA_BINARIES[config["system"]]
    except KeyError:
        sys.exit(f"❌ Unknown system '{config.get('system')}'. Must be one of: {list(SRA_BINARIES.keys())}")

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
            "snp_svant"
        shell:
            '''
            mkdir -p {params.out_dir}
            mkdir -p {params.data_dir}
            mkdir -p {config[logs]}
            echo {wildcards.sample}
            echo "Downloading {wildcards.sample} with {threads} threads"
            {FASTERQ} {wildcards.sample} -e {threads} 2>&1 {log}
            mv {wildcards.sample}*.fastq {params.data_dir}
            '''

