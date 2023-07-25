import os
import pandas as pd
from snakemake.utils import min_version
from snakemake.io import expand, glob_wildcards


##### set minimum snakemake version #####
min_version("6.10")


__author__ = ['Deepika Gunasekaran']
__version__ = '1.0.0'

configfile: "config/config.yaml"
report: "workflow/report/workflow.rst"

##### Include Rules #####

include: "workflow/rules/data_dump.smk"
include: "workflow/rules/trimming.smk"
include: "workflow/rules/genome_build.smk"
include: "workflow/rules/align.smk"
include: "workflow/rules/samtools_sort.smk"
include: "workflow/rules/common.smk"

##### Import samples based on config file #####

input_fp = config["samples"]["metadata"]
pd_header = 0 if config["samples"]["header"] else None

def get_samples():
    input_pd = pd.read_csv(input_fp,sep="\t",header=pd_header,names=["sample_name"])
    return input_pd

# Set SAMPLES based on config file
if not config["import"]:
    wc = glob_wildcards(os.path.join(config["outdir"], "raw/{sample}_1.fastq"))
    SAMPLES = [item for sublist in wc for item in sublist]
else:
    SAMPLES = get_samples()["sample_name"].tolist()

INDEX = list(range(1,5))

##### Target Rules #####

if config["import"] and config["trimming"]["trim"]:
    rule all:
        input:
            expand(
                os.path.join(config["outdir"], "raw/{sample}_{read}.fastq"),
                sample=SAMPLES,
                read=[1, 2]
            ),
            expand(os.path.join(config["outdir"], "preprocessed/trimmomatic/{sample}_trimmed_{read}.fastq"),
                sample=SAMPLES,
                read=[1, 2]
            ),
            expand(
                os.path.join(config["outdir"],"preprocessed/mapped/{sample}_aligned.bam"),
                sample=SAMPLES
            ),
            expand(
                os.path.join(config["logs"],"bowtie2_metrics_{sample}.txt"),
                sample=SAMPLES
            ),
            expand(
                os.path.join(config["outdir"],"preprocessed/mapped/{sample}_aligned_sorted.bam"),
                sample=SAMPLES
            ),
            expand(os.path.join(os.path.dirname(config["reference"]["genome"]),
                os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + ".{index}.bt2"),
                index=INDEX
            ),
            expand(os.path.join(os.path.dirname(config["reference"]["genome"]),
                os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + ".rev.{rev_index}.bt2"),
                rev_index=[1, 2]
            )
elif not config["import"] and config["trimming"]["trim"]:
    rule all:
        input:
            expand(os.path.join(config["outdir"], "preprocessed/trimmomatic/{sample}_trimmed_{read}.fastq"),
                sample=SAMPLES,
                read=[1, 2]
            ),
            expand(
                os.path.join(config["outdir"],"preprocessed/mapped/{sample}_aligned.bam"),
                sample=SAMPLES
            ),
            expand(
                os.path.join(config["logs"],"bowtie2_metrics_{sample}.txt"),
                sample=SAMPLES
            ),
            expand(
                os.path.join(config["outdir"],"preprocessed/mapped/{sample}_aligned_sorted.bam"),
                sample=SAMPLES
            ),
            expand(os.path.join(os.path.dirname(config["reference"]["genome"]),
                os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + ".{index}.bt2"),
                index=INDEX
            ),
            expand(os.path.join(os.path.dirname(config["reference"]["genome"]),
                os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + ".rev.{rev_index}.bt2"),
                rev_index=[1, 2]
            )
elif config["import"] and not config["trimming"]["trim"]:
    rule all:
        input:
            expand(
                os.path.join(config["outdir"], "raw/{sample}_{read}.fastq"),
                sample=SAMPLES,
                read=[1, 2]
            ),
            expand(
                os.path.join(config["outdir"],"preprocessed/mapped/{sample}_aligned.bam"),
                sample=SAMPLES
            ),
            expand(
                os.path.join(config["logs"],"bowtie2_metrics_{sample}.txt"),
                sample=SAMPLES
            ),
            expand(
                os.path.join(config["outdir"],"preprocessed/mapped/{sample}_aligned_sorted.bam"),
                sample=SAMPLES
            ),
            expand(os.path.join(os.path.dirname(config["reference"]["genome"]),
                os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + ".{index}.bt2"),
                index=INDEX
            ),
            expand(os.path.join(os.path.dirname(config["reference"]["genome"]),
                os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + ".rev.{rev_index}.bt2"),
                rev_index=[1, 2]
            )
else:
    rule all:
        input:
            expand(
                os.path.join(config["input_dir"],"{sample}_{read}.fastq"),
                sample=SAMPLES,
                read=[1, 2]
            ),
            expand(
                os.path.join(config["outdir"],"preprocessed/mapped/{sample}_aligned.bam"),
                sample=SAMPLES
            ),
            expand(
                os.path.join(config["logs"],"bowtie2_metrics_{sample}.txt"),
                sample=SAMPLES
            ),
            expand(
                os.path.join(config["outdir"],"preprocessed/mapped/{sample}_aligned_sorted.bam"),
                sample=SAMPLES
            )

