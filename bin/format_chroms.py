#!/usr/bin/env python3
"""
Script to change the format of chromosome names in biological files
"""
import argparse
import fileinput
import re
import roman

def main(args):
    """Main script"""
    with fileinput.input(args.input) as infile:
        for line in infile:
            if line[0] == args.comment:
                print(line, end='')
            else:
                spl = line.strip().split()
                spl[args.column] = format_chrom(spl[args.column], args.chr, args.roman)
                print(*spl, sep=args.sep)


def format_chrom(chrom, prefix='chr', rome=False):
    """Reformat a chromosome name"""
    num = re.sub('[Cc]hromosome|[Cc]hr', '', chrom.strip())
    if rome:
        if not roman.is_roman(num):
            num = roman.int_to_roman(int(num))
    else:
        if roman.is_roman(num):
            num = roman.roman_to_int(num)

    return prefix + str(num)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input', metavar='I', type=str,
                        help="File in which to reformat chromosome names")

    parser.add_argument('--chr', '-c', type=str, default='chr',
                        help="Prefix to append to chromsome number")

    parser.add_argument('--roman', '-r', action='store_true',
                        help="Format chromosme numbers as Roman numerals")

    parser.add_argument('--sep', '-s', default=None,
                        help="Column separater")

    parser.add_argument('--comment', '-k', default='#',
                        help="Comment character, any lines starting with this are ignored")

    parser.add_argument('--column', '-l', default=0, type=int,
                        help="Column containing chromosome (base 0)")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
