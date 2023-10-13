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
    vcf_file = "../../test/SRR7801919_filtered_snps_final_r2.vcf"

    # Define the MQ threshold (you can adjust this value)
    mq_threshold = 40
    qd_threshold = 2
    fs_threshold = 60
    sor_threshold = 4

    # Initialize lists to store MQ values for passed and failed variants
    mq_passed = []
    mq_failed = []
    qd_passed = []
    qd_failed = []
    fs_passed = []
    fs_failed = []
    sor_passed = []
    sor_failed = []

    # Open the VCF file
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))

    # Iterate through the records in the VCF file
    for record in vcf_reader:
        # Check if the variant passes the MQ threshold
        if 'MQ' in record.INFO and record.INFO['MQ'] >= mq_threshold:
            mq_passed.append(record.INFO['MQ'])
        else:
            # Variants that do not pass the threshold
            mq_failed.append(record.INFO.get('MQ', 0))  # Default to 0 if MQ is missing

        if 'QD' in record.INFO and record.INFO['QD'] >= qd_threshold:
            qd_passed.append(record.INFO['QD'])
        else:
            # Variants that do not pass the threshold
            qd_failed.append(record.INFO.get('QD', 0))  # Default to 0 if MQ is missing

        if 'FS' in record.INFO and record.INFO['FS'] <= fs_threshold:
            fs_passed.append(record.INFO['FS'])
        else:
            # Variants that do not pass the threshold
            fs_failed.append(record.INFO.get('FS', 0))  # Default to 0 if MQ is missing

        if 'SOR' in record.INFO and record.INFO['SOR'] <= sor_threshold:
            sor_passed.append(record.INFO['SOR'])
        else:
            # Variants that do not pass the threshold
            sor_failed.append(record.INFO.get('SOR', 0))  # Default to 0 if MQ is missing


    # Create a single plot with different colors for passed and failed variants
    plt.figure(figsize=(10, 6))

    # Plot MQ values for passed variants in green
    plt.hist(mq_passed, bins=20, color='blue', alpha=0.7, label='Passed Variants')

    # Plot MQ values for failed variants in red
    plt.hist(mq_failed, bins=20, color='red', alpha=0.7, label='Failed Variants')

    plt.xlabel('Mapping Quality (MQ)')
    plt.ylabel('Frequency')
    plt.title('MQ Distribution for Passed and Failed Variants')
    plt.legend()

    # Display the plot
    #plt.show()
    plt.savefig('../../results/MQ_filter_round2.png', dpi=300, bbox_inches='tight')

    # Create a single plot with different colors for passed and failed variants
    plt.figure(figsize=(10, 6))

    # Plot QD values for passed variants in green
    plt.hist(qd_passed, bins=20, color='blue', alpha=0.7, label='Passed Variants')

    # Plot MQ values for failed variants in red
    plt.hist(qd_failed, bins=20, color='red', alpha=0.7, label='Failed Variants')

    plt.xlabel('Quality normalized by depth (QD)')
    plt.ylabel('Frequency')
    plt.title('QD Distribution for Passed and Failed Variants')
    plt.legend()

    # Display the plot
    # plt.show()
    plt.savefig('../../results/QD_filter_round2.png', dpi=300, bbox_inches='tight')

    # Create a single plot with different colors for passed and failed variants
    plt.figure(figsize=(10, 6))

    # Plot QD values for passed variants in green
    plt.hist(fs_passed, bins=20, color='blue', alpha=0.7, label='Passed Variants')

    # Plot MQ values for failed variants in red
    plt.hist(fs_failed, bins=20, color='red', alpha=0.7, label='Failed Variants')

    plt.xlabel('Strand bias (FS)')
    plt.ylabel('Frequency')
    plt.title('FS Distribution for Passed and Failed Variants')
    plt.legend()

    # Display the plot
    # plt.show()
    plt.savefig('../../results/FS_filter_round2.png', dpi=300, bbox_inches='tight')

    # Create a single plot with different colors for passed and failed variants
    plt.figure(figsize=(10, 6))

    # Plot QD values for passed variants in green
    plt.hist(sor_passed, bins=20, color='blue', alpha=0.7, label='Passed Variants')

    # Plot MQ values for failed variants in red
    plt.hist(sor_failed, bins=20, color='red', alpha=0.7, label='Failed Variants')

    plt.xlabel('Strand odds ratio (SOR)')
    plt.ylabel('Frequency')
    plt.title('SOR Distribution for Passed and Failed Variants')
    plt.legend()

    # Display the plot
    # plt.show()
    plt.savefig('../../results/SOR_filter_round2.png', dpi=300, bbox_inches='tight')

    return


if __name__ == "__main__":
    main()
