#!/usr/bin/env python3
"""
Extract strains carrying each mutation in a VCF file
"""
import argparse
import fileinput
import vcf
from format_chroms import format_chrom

def main(args):
    """Main script"""
    # Import strain info
    if args.meta:
        meta = {}
        with fileinput.input(args.meta) as meta_file:
            for line in meta_file:
                if meta_file.isfirstline():
                    header = line.strip().split('\t')
                    standard = header.index('Standardized name')
                    isolate = header.index('Isolate name')

                else:
                    line = line.strip().split('\t')
                    meta[line[standard]] = line[isolate]

    vcf_file = vcf.Reader(filename=args.vcf)

    # Print headers
    if args.meta:
        print('mut_id', *[meta[i] for i in vcf_file.samples], sep='\t')

    else:
        print('mut_id', *vcf_file.samples, sep='\t')

    num_chroms = 2 * len(vcf_file.samples)
    # Iterate over loci
    for site in vcf_file:
        for i in range(len(site.ALT)):
            mut_id = ''.join((format_chrom(site.CHROM, rome=True),
                              ':', str(site.POS), '_', site.REF,
                              '/', str(site.ALT[i])))

            gens = []
            for call in site.samples:
                gens.append(call.data.GT.count(str(i + 1)))

            # Test frequency of variant and print
            if not args.filter or sum(gens) / num_chroms < args.filter:
                print(mut_id, *gens, sep='\t')

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('vcf', metavar='V', help="Input VCF file")

    parser.add_argument('--meta', '-m', default='',
                        help="File giving strain details to convert strain names to \
                              standardized names used in VCF")

    parser.add_argument('--filter', '-f', default=0.0, type=float,
                        help="Filter to genotypes at frequencies lower than the given level")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
