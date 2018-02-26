#!/usr/bin/env python3
"""
Extract mutfunc ID list from a VCF file
"""
import argparse
import fileinput
from format_chroms import format_chrom

def main(args):
    """Main script"""
    print('mut_id')
    with fileinput.input(args.vcf) as vcf_file:
        for line in vcf_file:
            if not line[0] == '#':
                spl = line.strip().split()
                if ',' in spl[4]:
                    for i in spl[4].split(','):
                        print(format_chrom(spl[0], args.prefix, args.roman),
                              ':', spl[1], '_', spl[3], '/', i, sep='')
                else:
                    print(format_chrom(spl[0], args.prefix, args.roman),
                          ':', spl[1], '_', spl[3], '/', spl[4], sep='')

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('vcf', metavar='V', type=str,
                        help="VCF input file")

    parser.add_argument('--roman', '-r', action='store_true',
                        help="Output chromosome numbers as Roman numerals")

    parser.add_argument('--prefix', '-p', default='chr',
                        help="Chromosome prefix")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
