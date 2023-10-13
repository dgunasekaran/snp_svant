#!/usr/bin/env python3.9


__author__ = "Deepika Gunasekaran"
__version__ = "0.0.1"
__maintainer__ = "Deepika Gunasekaran"
__email__ = "dgunasekaran@ucmerced.edu"
__status__ = "Development"


import os
import vcf
import argparse
import pandas as pd
from plotnine import *
import seaborn as sns
import matplotlib.pyplot as plt


def extract_values_from_vcf(vcf_path, info_field):
    vcf_reader = vcf.Reader(open(vcf_path, 'r'))
    values = []
    for record in vcf_reader:
        if info_field in record.INFO:
            values.append(record.INFO[info_field])
    return values


def main():
    # Replace 'your_file.vcf' with the path to your VCF file
    #vcf_file = vcf.Reader(open('../../test/preprocessed/SRR7801919_filtered_snps_final.vcf', 'r'))
    vcf_file = "../../test/SRR7801919_raw_snps_recal.vcf"

    # Thresholds
    mq_threshold = 40
    qd_threshold = 2
    fs_threshold = 60
    sor_threshold = 4
    mqrs_threshold = -12.5
    rprs_threshold = -8

    # Extract the values for MQ, QD, and FS
    mq_values = extract_values_from_vcf(vcf_file, 'MQ')
    qd_values = extract_values_from_vcf(vcf_file, 'QD')
    fs_values = extract_values_from_vcf(vcf_file, 'FS')
    sor_values = extract_values_from_vcf(vcf_file, 'SOR')
    mqrs_values = extract_values_from_vcf(vcf_file, 'MQRankSum')
    rprs_values = extract_values_from_vcf(vcf_file, 'ReadPosRankSum')

    # Create subplots for the three distributions
    fig, axes = plt.subplots(2, 3, figsize=(12, 6))

    # Plot MQ distribution with a vertical line
    sns.histplot(mq_values, kde=True, ax=axes[0, 0])
    axes[0, 0].axvline(x=mq_threshold, color='red', linestyle='--',
                       label=f'MQ > {mq_threshold} are retained')
    axes[0, 0].set_xlabel('MQ Values')
    axes[0, 0].set_ylabel('Density')
    axes[0, 0].set_title('MQ Distribution')
    axes[0, 0].legend()

    # Plot QD distribution with a vertical line
    sns.histplot(qd_values, kde=True, ax=axes[0, 1])
    axes[0, 1].axvline(x=qd_threshold, color='red', linestyle='--',
                       label=f'QD > {qd_threshold} are retained')
    axes[0, 1].set_xlabel('QD Values')
    axes[0, 1].set_ylabel('Density')
    axes[0, 1].set_title('QD Distribution')
    axes[0, 1].legend()

    # Plot FS distribution with a vertical line
    sns.histplot(fs_values, kde=True, ax=axes[0, 2])
    axes[0, 2].axvline(x=fs_threshold, color='red', linestyle='--',
                       label=f'FS < {fs_threshold} are retained')
    axes[0, 2].set_xlabel('FS Values')
    axes[0, 2].set_ylabel('Density')
    axes[0, 2].set_title('FS Distribution')
    axes[0, 2].legend()

    # Plot SOR distribution with a vertical line
    sns.histplot(sor_values, kde=True, ax=axes[1, 0])
    axes[1, 0].axvline(x=sor_threshold, color='red', linestyle='--',
                       label=f'SOR < {sor_threshold} are retained')
    axes[1, 0].set_xlabel('SOR Values')
    axes[1, 0].set_ylabel('Density')
    axes[1, 0].set_title('SOR Distribution')
    axes[1, 0].legend()

    # Plot MQRankSum distribution with a vertical line
    sns.histplot(mqrs_values, kde=True, ax=axes[1, 1])
    axes[1, 1].axvline(x=mqrs_threshold, color='red', linestyle='--',
                       label=f'MQRankSum > {mqrs_threshold} are retained')
    axes[1, 1].set_xlabel('MQRankSum Values')
    axes[1, 1].set_ylabel('Density')
    axes[1, 1].set_title('MQRankSum Distribution')
    axes[1, 1].legend()

    # Plot ReadPosRankSum distribution with a vertical line
    sns.histplot(rprs_values, kde=True, ax=axes[1, 2])
    axes[1, 2].axvline(x=rprs_threshold, color='red', linestyle='--',
                       label=f'ReadPosRankSum > {rprs_threshold} are retained')
    axes[1, 2].set_xlabel('ReadPosRankSum Values')
    axes[1, 2].set_ylabel('Density')
    axes[1, 2].set_title('ReadPosRankSum Distribution')
    axes[1, 2].legend()

    # Adjust layout and display the plot
    plt.tight_layout()
    #plt.show()
    plt.savefig('../../results/Figure_3.png', dpi=300, bbox_inches='tight')




    return


if __name__ == "__main__":
    main()
