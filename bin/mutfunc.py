#!/usr/bin/env python3
"""
Link functionally annotated and genomic mutations with their precalculated consequences

ToDo
- Rework import to be more versatile
- fill in data from mut id when importing effect file
- genomic as default type?
"""
import argparse
import fileinput

def main(args):
    """Main script"""
    mutations = {}

    # Import functional mutations
    if args.effect:
        with fileinput.input(args.functional) as fun_file:
            header = next(fun_file).strip().split('\t')
            for line in fun_file:
                line = line.strip().split('\t')
                line_id = ' '.join(line[:2])
                mutations[line_id] = {'effect': line[2],
                                      'gene': line[0],
                                      'pos_aa': line[1][1:-1],
                                      'ref_aa': line[1][:1],
                                      'alt_aa': line[1][-1]}

    else:
        with fileinput.input(args.functional) as fun_file:
            header = next(fun_file).strip().split('\t')
            for line in fun_file:
                line = line.strip().split('\t')
                line_id = ''.join((line[6], ' ', line[4], line[3], line[5]))
                mutations[line_id] = {header[i]:line[i] for i in range(8)}

    # Process SIFT data
    add_to_dict(mutations, ''.join((args.mutfunc, 'sift.tsv')),
                10, {'acc':0, 'sift_score':4, 'sift_median_ic':5})

    # Process Elms data
    add_to_dict(mutations, ''.join((args.mutfunc, 'elms.tsv')),
                12, {'elm_lost':9})

    # Process foldx data
    add_to_dict(mutations, ''.join((args.mutfunc, 'foldx.tsv')),
                0, {'foldx_ddG':3, 'foldx_ddG_sd':4, 'foldx_evidence':5})

    # Process foldx interaction data
    add_to_dict(mutations, ''.join((args.mutfunc, 'foldx-int.tsv')),
                0, {'foldx_int_ddG':7, 'foldx_int_ddG_sd':8,
                    'foldx_int_evidence':6, 'foldx_int_interactor':3})

    # Process phosphorylation data
    add_to_dict(mutations, ''.join((args.mutfunc, 'pho.tsv')),
                17, {'pho_prob':15})

    # Process PTM data
    add_to_dict(mutations, ''.join((args.mutfunc, 'ptms.tsv')),
                11, {'ptm_modification':8})

    # Process Genomic mutations
    # bit of a hack to work with the effect file which has no genomic ids
    if not args.effect:
        genomic_to_id = {}
        for key, value in mutations.items():
            genomic_to_id[value['mut_id']] = key

        if args.genomic:
            genomic_mutations = {}
            with fileinput.input(args.genomic) as gen_file:
                next(gen_file)
                for line in gen_file:
                    line = line.strip()
                    if not line in genomic_to_id.keys():
                        genomic_mutations[line] = {'type': 'genomic'}

        # Process TF data
        with fileinput.input(''.join((args.mutfunc, 'tfbs.tsv'))) as tf_file:
            next(tf_file)
            for line in tf_file:
                line = line.strip().split()
                if line[29] in genomic_to_id.keys():
                    mutations[genomic_to_id[line[29]]]['tf_score_diff'] = line[28]
                    mutations[genomic_to_id[line[29]]]['tf_perc_diff'] = line[27]
                elif args.genomic and line[29] in genomic_mutations.keys():
                    genomic_mutations[line[29]]['tf_score_diff'] = line[28]
                    genomic_mutations[line[29]]['tf_perc_diff'] = line[27]


    # Print Results
    values = ('mut_id', 'ref_codon', 'alt_codon', 'pos_aa', 'ref_aa',
              'alt_aa', 'gene', 'type', 'acc', 'sift_score', 'sift_median_ic',
              'elm_lost', 'foldx_ddG', 'foldx_ddG_sd', 'foldx_evidence',
              'foldx_int_ddG', 'foldx_int_ddG_sd', 'foldx_int_evidence',
              'foldx_int_interactor', 'pho_prob', 'ptm_modification',
              'tf_score_diff', 'tf_perc_diff', 'effect')

    defaults = ('NA', 'NA', 'NA', 'NA', 'NA',
                'NA', 'NA', 'NA', 'NA', 'NA', 'NA',
                'False', 'NA', 'NA', 'NA',
                'NA', 'NA', 'NA',
                'NA', 'NA', 'NA',
                'NA', 'NA', 'NA')

    defaults = dict(zip(values, defaults))

    process_missing(mutations, defaults)

    print('id', *values, sep='\t')
    for key, value in mutations.items():
        print(key, *(value[k] for k in values), sep='\t')

    if args.genomic:
        process_missing(genomic_mutations, defaults)
        for key, value in genomic_mutations.items():
            print('NA', key, *(value[k] for k in values[1:]), sep='\t')

def add_to_dict(dictionary, path, id_field, fields, sep='\t'):
    """Add fields from a table file to a dictionary via an ID field"""
    with fileinput.input(path) as in_file:
        next(in_file)
        for line in in_file:
            line = line.strip().split(sep)
            if line[id_field] in dictionary.keys():
                for key, value in fields.items():
                    dictionary[line[id_field]][key] = line[value]

def process_missing(dictionary, defaults):
    """Add defaults to a dictionary."""
    for value in dictionary.values():
        for i in defaults.keys():
            if i not in value.keys():
                value[i] = defaults[i]

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('functional', metavar='F',
                        help="Functionally annotated mutations file")

    parser.add_argument('--genomic', '-g',
                        help="Genomic mutations file")

    parser.add_argument('--mutfunc', '-m', default='mutfunc/',
                        help="Compiled Mutfunc directory")

    parser.add_argument('--effect', '-e', action='store_true',
                        help="Functional file gives a list of mutations and effects")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
