#!/usr/bin/env python3
"""
Script to determine probability of hog network disruption based on gene impact probabilities
"""
import argparse
import fileinput
import pathway as path

def main(args):
    """Main script"""
    meta = {}
    with fileinput.input(args.genes) as meta_file:
        for pair in meta_file:
            pair = pair.strip().split('\t')
            meta[pair[1]] = pair[0].lower()

    # Define Hog network
    hog = path.ProteinNetwork()
    hog.node(path.ProteinComplex('sho1-hkr1-msb1',
                                 [path.Protein('sho1'), path.Protein('hkr1'),
                                  path.Protein('msb1')]
                                ))

    hog.node(path.Protein('sln1'))

    hog.node(path.ProteinComplex('cla4-ste20-cdc42',
                                 [path.Protein('cla4'), path.Protein('ste20'),
                                  path.Protein('cdc42')], inputs=['sho1-hkr1-msb1']
                                ))

    hog.node(path.ProteinComplex('ste11-ste50',
                                 [path.Protein('ste11'), path.Protein('ste50')],
                                 inputs=['cla4-ste20-cdc42']
                                ))

    hog.node(path.Protein('ypd1', inputs=['sln1']))
    hog.node(path.Protein('ssk1', inputs=['ypd1']))
    hog.node(path.ProteinComplex('ssk2-ssk22',
                                 [path.Protein('ssk2'), path.Protein('ssk22')],
                                 inputs=['ssk1']
                                ))

    hog.node(path.Protein('pbs2', inputs=['ste11-ste50', 'ssk2-ssk22']))
    hog.node(path.Protein('hog1', inputs=['pbs2']))
    # hog.node(Protein('hot1', inputs=['hog1']))
    # hog.node(Protein('smp1', inputs=['hog1']))
    # hog.node(Protein('sko1', inputs=['hog1']))
    # hog.node(Protein('msn2', inputs=['hog1']))
    # hog.node(Protein('msn4', inputs=['hog1']))

    hog.set_input('sho1-hkr1-msb1', 'sln1')
    hog.set_output('hog1') #'hot1', 'smp1', 'msn2', 'msn4')

    # Evaluate network in each strain
    with fileinput.input(args.probs) as prob_file:
        header = next(prob_file).strip().split('\t')

        # Determine which columns have hog pathway genes
        hog_cols = {meta[header[i]]:i for i in range(len(header))
                    if header[i] in meta.keys()
                    and meta[header[i]] in hog.get_protein_names()}

        print('strain', 'hog_active', 'hog_probability', sep='\t')
        for line in prob_file:
            line = line.strip().split('\t')
            # Set protein functions as P(!impact) then eval network
            for name, col in hog_cols.items():
                hog.set_probability(name, 1 - float(line[col]))

            hog.calc_activity_probabilities()
            path_activity = hog.get_activity()
            hog_prob = hog.nodes['hog1'].probability_active

            # Print strain, network function and p(impact) for other genes
            print(line[0], int(path_activity), hog_prob, sep='\t')


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
