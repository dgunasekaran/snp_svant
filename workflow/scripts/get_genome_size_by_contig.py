#!/usr/bin/env python2.7


__author__ = "Deepika Gunasekaran"
__version__ = "0.0.1"
__maintainer__ = "Deepika Gunasekaran"
__email__ = "dgunasekaran@ucmerced.edu"
__status__ = "Development"

# Title: Parse multi-sequence fasta file and split into individual sequences
# Description: This program uses as input a multi-sequence fasta file and outputs individual sequences as a fasta file
# in user-specified folder


import os
import argparse


def get_chrom_size(input_fp):
    input_fh = open(input_fp, 'r')
    input_lines = input_fh.readlines()
    contig_size_dict = {}
    for line in input_lines:
        if line.startswith('>'):
            chrom_name = line[1:].split()[0]
            contig_size_dict[chrom_name] = 0
        else:
            contig_size_dict[chrom_name] += len(line.strip())
    return contig_size_dict


def main():
    parser = argparse.ArgumentParser(prog="Get genome stats",
                                     description="Parse directory with genomes in fasta format and get genome size and "
                                                 "contig number")
    parser.add_argument('-i', '--input_fp', type=str,
                        help="Input fasta file path",
                        required=True)
    parser.add_argument('-o', '--output_fp', type=str,
                        help="Output text file path",
                        required=True)
    args = parser.parse_args()
    contig_dict = get_chrom_size(args.input_fp)
    output_fp = args.output_fp
    output_str = ""
    for contig in contig_dict.keys():
        output_str += "\t".join([contig, str(contig_dict[contig])]) + "\n"
    with open(output_fp, 'w') as output_fh:
        output_fh.write(output_str)
    return


if __name__ == "__main__":
    main()
