import os

if config["trimming"]["trim"]:
    rule trimming:
        input:
            read_pair_1 = os.path.join(config["input_dir"], "{sample}_1.fastq"),
            read_pair_2 = os.path.join(config["input_dir"], "{sample}_2.fastq")
        output:
            trimmed_1 = os.path.join(config["outdir"], "preprocessed/trimmomatic/{sample}_trimmed_1.fastq"),
            trimmed_2 = os.path.join(config["outdir"], "preprocessed/trimmomatic/{sample}_trimmed_2.fastq"),
            untrimmed_1=os.path.join(config["outdir"],"preprocessed/trimmomatic/{sample}_untrimmed_1.fastq"),
            untrimmed_2=os.path.join(config["outdir"],"preprocessed/trimmomatic/{sample}_untrimmed_2.fastq"),
        params:
            out_dir=config["outdir"],
            preprocess_dir=os.path.join(config["outdir"], "preprocessed/"),
            trimming_dir=os.path.join(config["outdir"], "preprocessed/trimmomatic/"),
            summary_file=os.path.join(config["logs"], "trimmomatic_{sample}_summary.txt"),
            adapters=config["trimming"]["adapters"],
            seed_mismatches=config["trimming"]["seedMismatches"],
            pe_clip=config["trimming"]["palindromeClipThreshold"],
            se_clip=config["trimming"]["simpleClipThreshold"],
            leading=config["trimming"]["leading"],
            trailing=config["trimming"]["trailing"],
            window_size=config["trimming"]["windowSize"],
            qual_threshold=config["trimming"]["requiredQuality"],
            min_length=config["trimming"]["minlength"]
        log:
            os.path.join(config["logs"], "trimmomatic_{sample}.log")
        threads:
            config["threads"]
        shell:
            '''
            mkdir -p {params.preprocess_dir}
            mkdir -p {params.trimming_dir}
            echo {wildcards.sample}
            trimmomatic PE -threads {threads} -phred33 -trimlog {log} -summary {params.summary_file} \
            {input.read_pair_1} {input.read_pair_2} \
            {output.trimmed_1} {output.untrimmed_1} {output.trimmed_2} {output.untrimmed_2} \
            ILLUMINACLIP:{params.adapters}:{params.seed_mismatches}:{params.pe_clip}:{params.se_clip} \
            LEADING:{params.leading} TRAILING:{params.trailing} \
            SLIDINGWINDOW:{params.window_size}:{params.qual_threshold} MINLEN:{params.min_length}
        '''
