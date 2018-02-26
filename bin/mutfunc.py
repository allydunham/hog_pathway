#!/usr/bin/env python3
"""
Link functionally annotated and genomic mutations with their precalculated consequences

ToDo
- Genomic Mutations
"""
import argparse
import fileinput

def main(args):
    """Main script"""
    mutations = {}

    # Import functional mutations
    with fileinput.input(args.functional) as fun_file:
        header = next(fun_file).strip().split('\t')
        for line in fun_file:
            line = line.strip().split('\t')
            line_id = ''.join((line[6], ' ', line[4], line[3], line[5]))
            mutations[line_id] = {header[i]:line[i] for i in range(8)}

    add_to_dict(mutations, ''.join((args.mutfunc, 'sift.tsv')),
                10, {'acc':0, 'sift_score':4, 'sift_median_ic':5})

    # # Process SIFT data
    # with fileinput.input(''.join((args.mutfunc, 'sift.tsv'))) as sift_file:
    #     next(sift_file)
    #     for line in sift_file:
    #         line = line.strip().split('\t')
    #         if line[10] in mutations.keys():
    #             mutations[line[10]]['acc'] = line[0]
    #             mutations[line[10]]['sift_score'] = line[4]
    #             mutations[line[10]]['sift_median_ic'] = line[5]

    # # Process Elms data
    # with fileinput.input(''.join((args.mutfunc, 'elms.tsv'))) as elms_file:
    #     next(elms_file)
    #     for line in elms_file:
    #         line = line.strip().split('\t')
    #         if line[12] in mutations.keys():
    #             mutations[line[12]]['elm_lost'] = line[9]

    # # Process foldx data
    # with fileinput.input(''.join((args.mutfunc, 'foldx.tsv'))) as foldx_file:
    #     next(foldx_file)
    #     for line in foldx_file:
    #         line = line.strip().split('\t')
    #         if line[0] in mutations.keys():
    #             mutations[line[0]]['foldx_ddG'] = line[3]
    #             mutations[line[0]]['foldx_ddG_sd'] = line[4]
    #             mutations[line[0]]['foldx_evidence'] = line[5]

    # # Process foldx interaction data
    # with fileinput.input(''.join((args.mutfunc, 'foldx-int.tsv'))) as foldx_int_file:
    #     next(foldx_int_file)
    #     for line in foldx_int_file:
    #         line = line.strip().split('\t')
    #         if line[0] in mutations.keys():
    #             mutations[line[0]]['foldx_int_ddG'] = line[7]
    #             mutations[line[0]]['foldx_int_ddG_sd'] = line[8]
    #             mutations[line[0]]['foldx_int_evidence'] = line[6]
    #             mutations[line[0]]['foldx_int_interactor'] = line[3]

    # # Process phosphorylation data
    # with fileinput.input(''.join((args.mutfunc, 'pho.tsv'))) as pho_file:
    #     next(pho_file)
    #     for line in pho_file:
    #         line = line.strip().split('\t')
    #         if line[17] in mutations.keys():
    #             mutations[line[17]]['pho_prob'] = line[15]

    # # Process PTM data
    # with fileinput.input(''.join((args.mutfunc, 'ptms.tsv'))) as ptm_file:
    #     next(ptm_file)
    #     for line in ptm_file:
    #         line = line.strip().split('\t')
    #         if line[11] in mutations.keys():
    #             mutations[line[11]]['ptm_modification'] = line[8]

    # Process Genomic mutations
    genomic_to_id = {}
    for key, value in mutations.items():
        genomic_to_id[value['mut_id']] = key

    genomic_mutations = {}
    with fileinput.input(args.genomic) as gen_file:
        next(gen_file)
        for line in gen_file:
            line = line.strip()
            if not line in genomic_to_id.keys():
                genomic_mutations[line] = {}

    # Process TF data
    with fileinput.input(''.join((args.mutfunc, 'tfbs.tsv'))) as tf_file:
        next(tf_file)
        for line in tf_file:
            line = line.strip().split()
            if line[29] in genomic_to_id.keys():
                mutations[genomic_to_id[line[29]]]['tf_score_diff'] = line[28]
                mutations[genomic_to_id[line[29]]]['tf_perc_diff'] = line[27]
            elif line[29] in genomic_mutations.keys():
                genomic_mutations[line[29]]['tf_score_diff'] = line[28]
                genomic_mutations[line[29]]['tf_perc_diff'] = line[27]


    # Print Results
    vals = ('mut_id', 'ref_codon', 'alt_codon', 'pos_aa', 'ref_aa',
            'alt_aa', 'gene', 'type', 'acc', 'sift_score', 'sift_median_ic',
            'elm_lost', 'foldx_ddG', 'foldx_ddG_sd', 'foldx_evidence',
            'foldx_int_ddG', 'foldx_int_ddG_sd', 'foldx_int_evidence',
            'foldx_int_interactor', 'pho_prob', 'ptm_modification',
            'tf_score_diff', 'tf_perc_diff')

    print('id', *vals, sep='\t')
    for key, value in mutations.items():
        for i in vals:
            if not i in value.keys():
                value[i] = 'NA'

        print(key, *(value[k] for k in vals), sep='\t')

    for key, value in genomic_mutations.items():
        for i in vals[1:]:
            if not i in value.keys():
                value[i] = 'NA'

        print('NA', key, *(value[k] for k in vals[1:]), sep='\t')

def add_to_dict(dictionary, path, id_field, fields, sep='\t'):
    """Add fields from a table file to a dictionary via an ID field"""
    with fileinput.input(path) as in_file:
        next(in_file)
        for line in in_file:
            line = line.strip().split(sep)
            if line[id_field] in dictionary.keys():
                for key, value in fields.items():
                    dictionary[line[id_field]][key] = line[value]

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('functional', metavar='F',
                        help="Functionally annotated mutations file")

    parser.add_argument('genomic', metavar='G',
                        help="Genomic mutations file")

    parser.add_argument('--mutfunc', '-m', default='mutfunc/',
                        help="Compiled Mutfunc directory")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
