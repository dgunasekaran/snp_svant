# ABOUT:

SNP-SVant is a variant calling workflow that can be used to identify SNPs and SVs in non-benchmarked organisms. The SNP-SVant workflow requires as input files paired-end short-read FASTQ files or SRA IDs listed in test_samples.tsv (example: test/metadata/test_samples.tsv), reference FASTA file and reference annotation GFF file. The outputs of SNP-SVant include a Variant Call Format (VCF) file with SNPs and small indels obtained using GATK, a VCF file with structural variants obtained using GRIDSS, and annotated SVs listed in BED format. Additionally, the scores used to assess the quality of the variant calls can be visualized and used to adjust the variant filtering parameters in GATK.

# DOWNLOAD REPOSITORY AND SETUP:

Download git repository using the following commands:

```
git clone https://github.com/dgunasekaran/snp_svant
```

Setup Conda environment and install R packages using the following commands:

```
cd snp_svant
conda env create -f workflow/envs/environment.yaml
conda activate snp_svant
Rscript workflow/scripts/install_r_packages.R
```

# USAGE:

Before performing variant calling, build genome indices using the command:

```
snakemake --cores <number_of_threads> build --use-conda
```

Modify parameters in the `config/config.yaml` file to refer to input and output paths, modify variant filteration criteria and refering to genome fasta and annotation files. The description of the parameters are given in `workflow/schemas/config.schemas.yaml`. Run the workflow using command:

```
snakemake --cores <number_of_threads> all --use-conda
```

To generate graphs of QC metrics from GATK, create conda environment and use the commands:

```
conda create -n vcf_qc python pyvcf pandas argparse seaborn
conda activate vcf_qc
python workflow/scripts/variant_quality_assessment.py -i <input_vcf> -o <output_vcf> 
```

# RESOURCES:

1. GATK best practices can be used to understand quality thresholds and can be found at https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows.

1. GATK is used for SNP and INDEL calling and details about this method can be found in the [original publication of GATK](https://pubmed.ncbi.nlm.nih.gov/20644199/) and publication of [best practices](https://pubmed.ncbi.nlm.nih.gov/25431634/).

1. GRIDSS is used for structural variant calling and details about this method can be found in the original publications of [GRIDSS](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5741059/) and [GRIDSS2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02423-x).


# TROUBLESHOOTING:

1. If ImportError occurs while importing required packages, update Mamba in the Conda environment or re-install Mamba using the following command:
```
conda install mamba -c conda-forge
```

2. If you encounter dependecy issues with R packages, you can use the locked conda environment (after installing conda-lock) as follows:
```
pip install conda-lock
conda-lock install --name snp_svant_lock workflow/envs/snp_svant-lock.yml
conda activate snp_svant_lock
```

3. If output files are incomplete or if run crashes unexpectedly, re-run the pipeline or use additional parameter `--rerun-incomplete` to regenerate incomplete files.

4. If output from VEP is not generated and SNPs are not annotated, it is likely due to missing PERL dependencies required to run VEP. You can refer to VEP installation [here](https://useast.ensembl.org/info/docs/tools/vep/index.html), to install PERL module dependencies or re-run VEP outside the workflow. 


