import os

rule gridss_build:
    input:
        reference=config["reference"]["genome_gridss"],
    output:
        idx=multiext(config["reference"]["genome_gridss"], ".amb", ".ann", ".bwt", ".dict", ".fai", ".pac", ".sa")
    params:
        extra="--jvmheap 1g"
    log:
        "log/gridss/setupreference.log"
    wrapper:
        "v2.2.1/bio/gridss/setupreference"
