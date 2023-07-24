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

include: "workflow/rules/data_dump_external.smk"
include: "workflow/rules/common.smk"

##### Import samples based on config file #####

input_fp = config["samples"]["metadata"]
pd_header = 0 if config["samples"]["header"] else None

def get_samples():
    input_pd = pd.read_csv(input_fp,sep="\t",header=pd_header,names=["sample_name"])
    return input_pd

# Set SAMPLES based on config file
if not config["import"]:
    SAMPLES = glob_wildcards(os.path.join(config["outdir"], "raw/{sample}_1.fastq"))
else:
    SAMPLES = get_samples()["sample_name"].tolist()

##### Target Rules #####

rule all:
    input:
        expand(
            os.path.join(config["outdir"], "raw/{sample}_{read}.fastq"),
            sample=SAMPLES,
            read=[1, 2]
        )
