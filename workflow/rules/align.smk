import os

rule reference_align:
    input:
        sample=[
            os.path.join(config["outdir"],"preprocessed/trimmomatic/{sample}_trimmed_1.fastq"),
            os.path.join(config["outdir"],"preprocessed/trimmomatic/{sample}_trimmed_2.fastq")
            ] if config["trimming"]["trim"] else [
            os.path.join(config["outdir"],"raw/{sample}_1.fastq"),
            os.path.join(config["outdir"],"raw/{sample}_2.fastq")
            ] if config["import"] else [
            os.path.join(config["input_dir"],"{sample}_1.fastq"),
            os.path.join(config["input_dir"],"{sample}_2.fastq")
            ],
        idx=multiext(
            os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0]),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2"
        ),
    output:
        os.path.join(config["outdir"], "preprocessed/mapped/{sample}_aligned.bam"),
        metrics=os.path.join(config["logs"], "bowtie2_metrics_{sample}.txt"),
    log:
        os.path.join(config["logs"], "bowtie2_align_{sample}.log"),
    params:
        extra="--local",
    threads:
        config["threads"]
    wrapper:
        "v2.2.1/bio/bowtie2/align"
# bowtie2 -q -p 32 --local -x /phd_16tb/calbicans_biofilms/resources/bowtie2_index/calbicans_21/C_albicans_SC5314_A21_current_chromosomes_without_mt -1 /phd_16tb/calbicans_biofilms/data/raw/${strain}_1.fastq -2 /phd_16tb/calbicans_biofilms/data/raw/${strain}_2.fastq -S /phd_16tb/calbicans_biofilms/data/preprocessed/bowtie2/${strain}_aligned.sam; done 2>&1 >> /phd_16tb/calbicans_biofilms/data/preprocessed/bowtie2_local.log
