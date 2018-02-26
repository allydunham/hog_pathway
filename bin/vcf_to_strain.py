#!/usr/bin/env python3
"""
Extract strains carrying each mutation in a VCF file

ToDo
- deal with multiple alt alleles
"""
import argparse
#import fileinput
import vcf
from format_chroms import format_chrom

def main(args):
    """Main script"""
    vcf_file = vcf.Reader(filename=args.vcf)
    genotypes = {}
    for site in vcf_file:
        for i in range(len(site.ALT)):
            mut_id = ''.join((format_chrom(site.CHROM, rome=True),
                              ':', str(site.POS), '_', site.REF,
                              '/', str(site.ALT[i])))

            genotypes[mut_id] = []
            for call in site.samples:
                genotypes[mut_id].append(call.data.GT.count(str(i + 1)))

    # Print matrix of strain genotypes
    print('mut_id', *vcf_file.samples, sep='\t')
    for mut_id, gens in genotypes.items():
        print(mut_id, *gens, sep='\t')

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('vcf', metavar='V', help="Input VCF file")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
