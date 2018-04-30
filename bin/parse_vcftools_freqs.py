#!/usr/bin/env python3
"""
Script to convert a VCFtools .frq file into a list of variant/frequency pairs
with mutfunc format variant identifier.

ToDo
- more info about different allele frequencies
"""
import sys
import argparse
import fileinput
from vcf_to_strain import format_chrom

def main(args):
    """Main script"""
    print('mut_id', 'freq', sep='\t', file=sys.stdout)
    with fileinput.input(args.freqs) as frq_file:
        next(frq_file)
        for line in frq_file:
            line = line.strip().split('\t')
            chrom = format_chrom(line[0], prefix='chr', rome=True)
            alleles = [x.split(':') for x in line[4:]] # Gives [[allele, freq], ...] pairs
            base_id = ''.join([chrom, ':', line[1], '_', alleles[0][0], '/'])
            for allele in alleles[1:]:
                print(base_id + allele[0], allele[1], sep='\t', file=sys.stdout)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('freqs', metavar='F', help="VCFtools .frq file")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
