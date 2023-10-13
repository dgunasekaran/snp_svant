#!/usr/bin/env python3.9


__author__ = "Deepika Gunasekaran"
__version__ = "0.0.1"
__maintainer__ = "Deepika Gunasekaran"
__email__ = "dgunasekaran@ucmerced.edu"
__status__ = "Development"


import os
import vcf
import gffutils
import tempfile
import argparse
from Bio import SeqIO


class Parser:
    @staticmethod
    def parse_arguments():
        parser = argparse.ArgumentParser(description="This program parses VCF files and returns aligned fasta")
        parser.add_argument('-i', '--input', type=str,
                            help='Input VCF file',
                            required=True)
        parser.add_argument('-r', '--reference', type=str, help="Reference fasta file used to generate VCF",
                            required=True)
        parser.add_argument('--include_ref', type=bool, help='Include reference in final aligned Fasta',
                            default=True,
                            required=False)
        parser.add_argument('--no_gaps', type=bool, help='Exclude gaps in final aligned fasta i.e. only positions '
                                                         'without any gaps in any of the samples will be included',
                            default=True,
                            required=False)
        parser.add_argument('-g', '--gff', type=str, help="GFF file to extract specific features for multi-fasta",
                            default=None,
                            required=False)
        parser.add_argument('-f', '--feature', type=str, help="Feature type to extract from GFF. Default: gene",
                            default='gene',
                            required=False)
        parser.add_argument('-l', '--locus', type=valid_locus, help='Extract fasta for specific locus. Should be of '
                                                                    'type chrom:start:end. The start and end positions '
                                                                    'are inclusive (1-indexed).',
                            required=False,
                            default=None)
        parser.add_argument('-s', '--strand', type=valid_strand, help='Strand should be + or -. Default: +',
                            required=False,
                            default='+')
        parser.add_argument('-o', '--output', type=str, help='Output filename or directory (if GFF provided with a list'
                                                             ' of features for multiple multi-fasta files).',
                            required=False)
        args = parser.parse_args()
        return args


def valid_locus(locus: str):
    try:
        chrom, star, end = str(locus.split(':')[0]), int(locus.split(':')[1]), int(locus.split(':')[2])
        return locus
    except (IndexError, ValueError, TypeError) as error:
        print(f"ERROR: {error}")
        raise argparse.ArgumentTypeError(f"Given locus: {locus} is not valid. Specify locus in format chrom:start:end")


def valid_strand(strand: str):
    if strand == '+' or strand == '-':
        return strand
    else:
        raise argparse.ArgumentTypeError(f"Given strand: {strand} is not valid. Strand should be '+' or '-'")


def get_reference_sequences(ref_fp):
    ref_dict = {}
    for record in SeqIO.parse(ref_fp, 'fasta'):
        ref_dict[record.id.split(' ')[0]] = str(record.seq)
    return ref_dict


def get_record_alleles(record, samples):
    sample_allele_frequency = {}
    record_alleles = [record.REF]
    for alt_allele in record.ALT:
        record_alleles.append(alt_allele)
    # Filter alleles to remove gaps i.e. insertion or deletion alternatives
    filtered_alleles = []
    filtered_allele_indices = []
    for ind in range(0, len(record_alleles)):
        allele = record_alleles[ind]
        if allele != '*' and len(allele) == 1:
            filtered_alleles.append(allele)
            filtered_allele_indices.append(ind)
    filtered_alleles_dict = dict(zip(filtered_allele_indices, filtered_alleles))
    for sample in samples:
        sample_name = sample.sample
        sample_call = record.genotype(sample_name)
        sample_allele_frequency[sample_name] = {}
        # If GT is ./. which denotes missing information, assume no variant
        for ind, allele in filtered_alleles_dict.items():
            if sample_call['AD'] is not None:
                sample_allele_frequency[sample_name][str(allele)] = sample_call['AD'][ind]
    return sample_allele_frequency


def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_sequence = sequence[::-1]
    reverse_complement_sequence = ''.join(complement_dict.get(base, base) for base in reverse_sequence)
    return reverse_complement_sequence


def get_alignments(vcf_fp, ref_seqs: dict, strand, locus=None):
    added_chrom = []
    aln_str = ""
    vcf_reader = vcf.Reader(filename=vcf_fp)
    vcf_records = []
    sample_chrom_dict = {}
    for record in vcf_reader:
        if locus is not None:
            l_chrom, l_start, l_end = str(locus.split(':')[0]), int(locus.split(':')[1]), int(locus.split(':')[2])
            if record.CHROM == l_chrom and (int(record.POS) >= l_start) and (int(record.POS) <= l_end):
                vcf_records.append(record)
            elif record.CHROM == l_chrom and int(record.POS) > l_end:
                break
        else:
            vcf_records.append(record)
    for rec in vcf_records:
        chrom = rec.CHROM
        if chrom not in ref_seqs.keys():
            print(f"Reference does not match VCF (chromosome {chrom} not found in reference). Ensure the reference "
                  f"fasta is the same as the one used for generating the vcf. Skipping variant!")
        else:
            reference_seq = ref_seqs[chrom]
            start_index = 1
            if locus is not None:
                ref_substring = reference_seq[l_start-1:l_end]
                reference_seq = ref_substring
                start_index = l_start
            if chrom not in added_chrom:
                chrom_header = f">Reference_{chrom}\n"
                if strand == '+':
                    aln_str += chrom_header + reference_seq + "\n"
                    added_chrom.append(chrom)
                else:
                    aln_str += chrom_header + reverse_complement(reference_seq) + "\n"
                    added_chrom.append(chrom)
            sample_af = get_record_alleles(rec, rec.samples)
            for sample in rec.samples:
                sample_name = sample.sample
                sample_chrom_seq = reference_seq
                if sample_name not in sample_chrom_dict.keys():
                    sample_chrom_dict[sample_name] = {chrom: reference_seq}
                elif chrom not in sample_chrom_dict[sample_name].keys():
                    sample_chrom_dict[sample_name][chrom] = reference_seq
                else:
                    sample_chrom_seq = sample_chrom_dict[sample_name][chrom]
                # If alleles in sample is not empty i.e. if alternate alleles are present for sample get allele with
                # maximum depth
                # LIMITATION: In case of heterozygous sites, only the allele with maximum depth with be included
                if sample_af[sample_name] != {}:
                    sample_af[sample_name] = {key: value if value is not None else 0
                                              for key, value in sample_af[sample_name].items()}
                    sample_max_allele = max(sample_af[sample_name], key=sample_af[sample_name].get)
                    ref_position = rec.POS
                    change_position = ref_position - start_index
                    sample_chrom_seq_w_variant = sample_chrom_seq[:change_position] + sample_max_allele + \
                                                 sample_chrom_seq[change_position+1:]
                    sample_chrom_dict[sample_name][chrom] = sample_chrom_seq_w_variant
    if strand != '+':
        rev_comp_dict = {}
        for sample in sample_chrom_dict.keys():
            rev_comp_dict[sample] = {}
            for chrom in sample_chrom_dict[sample].keys():
                rev_comp_dict[sample][chrom] = reverse_complement(sample_chrom_dict[sample][chrom])
        return aln_str, rev_comp_dict
    return aln_str, sample_chrom_dict


def get_output_fp(arg_output, arg_input, out_fn="vcf_to_aln.fasta"):
    if arg_output is None:
        print(f"Output directory not specified! Writing file to {os.path.join(os.path.dirname(arg_input), out_fn)}")
        return os.path.join(os.path.dirname(arg_input), out_fn)
    elif os.path.isdir(arg_output):
        return os.path.join(arg_output, out_fn)
    elif os.path.exists(arg_output):
        output_dir = os.path.dirname(arg_output)
        # Create a temporary file
        temp_file = tempfile.NamedTemporaryFile(dir=output_dir)
        # Get the path of the temporary file
        temp_file_path = temp_file.name
        print(f"File {arg_output} exists! Writing output to {temp_file_path}")
        return temp_file_path
    return arg_output


def write_fasta(header, sample_dict, output_fp, include_ref):
    if not include_ref:
        header = ""
    for sample in sample_dict.keys():
        for chrom in sample_dict[sample]:
            header += f">{sample}_{chrom}\n{sample_dict[sample][chrom]}\n"
    with open(output_fp, 'w') as output_fh:
        output_fh.write(header)
    return


def main():
    prog_args = Parser.parse_arguments()
    ref_seqs = get_reference_sequences(prog_args.reference)
    if prog_args.gff is None or (not os.path.exists(prog_args.gff)):
        if prog_args.gff is not None and (not os.path.exists(prog_args.gff)):
            print(f"Path for GFF {prog_args.gff} does not exist! Skipping GFF!")
        header, sample_seq_dict = get_alignments(vcf_fp=prog_args.input, ref_seqs=ref_seqs, strand=prog_args.strand,
                                                 locus=prog_args.locus)
        output_fp = get_output_fp(prog_args.output, prog_args.input)
        write_fasta(header, sample_seq_dict, output_fp, prog_args.include_ref)
    else:
        # Check for GFF
        gff_lines = open(prog_args.gff, 'r').readlines()
        all_line_elements = []
        for line in gff_lines:
            if line.startswith('#'):
                continue
            line_elements = line.strip().split('\t')
            for element in line_elements:
                all_line_elements.append(element)
        if prog_args.feature not in all_line_elements:
            print(f"The feature type '{prog_args.feature}' does not exist in the GFF file. Using 'gene' instead!")
        elif 'gene' not in all_line_elements:
            print(f"'gene' does not exist in the GFF file! GFF file could possibly be of incorrect format. Please "
                  f"specify correct feature to use or ensure GFF file contains 'gene' specification! Skiping GFF!")
            header, sample_seq_dict = get_alignments(vcf_fp=prog_args.input, ref_seqs=ref_seqs, strand=prog_args.strand,
                                                     locus=prog_args.locus)
            output_fp = get_output_fp(prog_args.output, prog_args.input)
            write_fasta(header, sample_seq_dict, output_fp, prog_args.include_ref)
        else:
            db = gffutils.create_db(prog_args.gff, dbfn=':memory:', force=True, id_spec=prog_args.feature)
            if prog_args.locus is not None:
                l_chrom, l_start, l_end = str(prog_args.locus.split(':')[0]), int(prog_args.locus.split(':')[1]), \
                    int(prog_args.locus.split(':')[2])
            # Iterate over the features in the GFF file
            for feature in db.features_of_type(prog_args.feature):
                # Access feature attributes
                feature_name = feature.attributes['Name'][0]
                if not feature_name.startswith('orf'):
                    continue
                feature_start = feature.start
                feature_end = feature.end
                feature_chrom = feature.seqid
                feature_strand = feature.strand
                if feature_chrom != l_chrom:
                    continue
                elif feature_start >= l_start and feature_end <= l_end:
                    locus = ':'.join([feature_chrom, str(feature_start), str(feature_end)])
                    output_fn = feature_name + '.fasta'
                    header, sample_seq_dict = get_alignments(vcf_fp=prog_args.input, ref_seqs=ref_seqs,
                                                             strand=feature_strand,
                                                             locus=locus)
                    output_fp = get_output_fp(prog_args.output, prog_args.input, out_fn=output_fn)
                    write_fasta(header, sample_seq_dict, output_fp, prog_args.include_ref)
        # TO DO: Include GFF and locus specification i.e. extract features only from specific loci
    return


if __name__ == "__main__":
    main()
