import os

rule bowtie_build:
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
    conda:
        "../envs/bowtie2.yaml"
    wrapper:
        "v2.2.1/bio/bowtie2/build"
