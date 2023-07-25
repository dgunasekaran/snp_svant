import os

rule bowtie2_build:
    input:
        ref=os.path.abspath(config["reference"]["genome"])
    output:
        os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.1.bt2'),
        os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.2.bt2'),
        os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.3.bt2'),
        os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.4.bt2'),
        os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.rev.1.bt2'),
        os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.rev.2.bt2')
    log:
        os.path.join(config["logs"], "bowtie2_build.log")
    threads:
        config["threads"]
    wrapper:
        "v2.2.1/bio/bowtie2/build"

rule fasta_index:
    input:
        os.path.abspath(config["reference"]["genome"]),
    output:
        os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.fasta.fai'),
    log:
        os.path.join(config["logs"],"samtools_fasta_index.log")
    wrapper:
        "v2.2.1/bio/samtools/faidx"


rule genome_size:
    input:
        os.path.abspath(config["reference"]["genome"]),
    output:
        os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.txt'),
    shell:
        '''
        python workflow/scripts/get_genome_size_by_contig.py -i {input} -o {output}
        '''


rule picard_create_dict:
    input:
        os.path.abspath(config["reference"]["genome"]),
    output:
        os.path.join(os.path.dirname(config["reference"]["genome"]),
            os.path.splitext(os.path.basename(config["reference"]["genome"]))[0] + '.dict'),
    log:
        os.path.join(config["logs"],"picard_create_dict.log"),
    resources:
        mem_mb=1024,
    wrapper:
        "v2.2.1/bio/picard/createsequencedictionary"
