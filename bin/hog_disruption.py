#!/usr/bin/env python3
"""
Script to determine probability of hog network disruption based on gene impact probabilities
"""
import argparse
import fileinput
import binary_pathway as bp

def main(args):
    """Main script"""
    meta = {}
    with fileinput.input(args.genes) as meta_file:
        for pair in meta_file:
            pair = pair.strip().split('\t')
            meta[pair[1]] = pair[0].lower()

    # Define Hog network
    hog = bp.ProteinNetwork()
    hog.node(bp.Protein('cdc42'))
    hog.node(bp.ProteinComplex('sho1-hkr1-msb1',
                               [bp.Protein('sho1'), bp.Protein('hkr1'), bp.Protein('msb1')]))

    hog.node(bp.Protein('sln1'))
    hog.node(bp.ProteinComplex('cla4-ste20',
                               [bp.Protein('cla4'), bp.Protein('ste20')],
                               inputs=['cdc42']))

    hog.node(bp.Protein('ste11', inputs=['cla4-ste20', 'sho1-hkr1-msb1']))
    hog.node(bp.Protein('ypd1', inputs=['sln1']))
    hog.node(bp.Protein('ssk1', inputs=['ypd1']))
    hog.node(bp.ProteinComplex('ssk2-ssk22',
                               [bp.Protein('ssk2'), bp.Protein('ssk22')],
                               inputs=['ssk1']))

    hog.node(bp.Protein('pbs2', inputs=['ste11', 'ssk2-ssk22']))
    hog.node(bp.Protein('hog1', inputs=['pbs2']))
    # hog.node(Protein('hot1', inputs=['hog1']))
    # hog.node(Protein('smp1', inputs=['hog1']))
    # hog.node(Protein('sko1', inputs=['hog1']))
    # hog.node(Protein('msn2', inputs=['hog1']))
    # hog.node(Protein('msn4', inputs=['hog1']))

    hog.set_input('cdc42', 'sho1-hkr1-msb1', 'sln1')
    hog.set_output('hog1') #'hot1', 'smp1', 'msn2', 'msn4')

    # Evaluate network in each strain
    with fileinput.input(args.probs) as prob_file:
        header = next(prob_file).strip().split('\t')
        print([meta[header[i]] for i in range(len(header)) if header[i] in meta.keys()])

        # Determine which columns have hog pathway genes
        gene_cols = {meta[header[i]]:i for i in range(len(header))
                     if header[i] in meta.keys()
                     and meta[header[i]] in hog.get_protein_names()}
        
        print('strain',sep='\t')
        for line in prob_file:
            line = line.strip().split('\t')
            # Set protein functions as P(!impact) then eval network
            for name, col in gene_cols.items():
                hog.set_function(name, 1 - float(line[col]))
            path_activity = hog.get_activity()

            # Print strain, network function and p(impact) for other genes
            print(line[0], int(path_activity), sep='\t')


def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('probs', metavar='P', help="Table of impact probabilities")

    parser.add_argument('--genes', '-g', default='meta/hog-gene-ids',
                        help="Table of gene name/id pairs")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
