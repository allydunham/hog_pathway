#!/usr/bin/env python3
"""
Process per mutation genotypes to group them by gene and assign per gene disruption scores

ToDo
- Currently nothing takes account of heterozygotes
- further logic around nonsense/frameshift
"""
import argparse
import fileinput
import sys
import numpy as np
import pandas as pd

def main(args):
    """Main script"""
    impacts = pd.read_table(args.mutations)

    genes = impacts.gene.dropna().unique()

    with fileinput.input(args.genotypes) as geno_file:
        header = next(geno_file).strip().split('\t')

        # Read in and filter strains
        if args.strains:
            try:
                strains = []
                with fileinput.input(args.strains) as strain_file:
                    for line in strain_file:
                        strains.append(line.strip())

            except FileNotFoundError:
                strains = args.strains.strip().split(',')

        else:
            strains = header[1:]

        indeces = [i for i in range(len(header)) if header[i] in strains]
        strain_muts = {s:{g:[] for g in genes} for s in strains}

        # Extract per gene per strain mutations
        for line in geno_file:
            line = line.strip().split('\t')
            mut_id = line[0]
            # Identify impacted genes
            affected_genes = impacts[impacts['mut_id'] == mut_id].gene.dropna().values
            for i in indeces:
                if not line[i] == '0':
                    for gene in affected_genes:
                        strain_muts[header[i]][gene].append(mut_id)

    # Determine per gene mutation probabilities
    print("strain", *genes, sep='\t')
    for strain, genes in strain_muts.items():
        probs = {}
        for gene, muts in genes.items():
            probs[gene] = gene_impact_prob(muts, impacts, gene)

        print(strain, *[probs[i] for i in genes], sep='\t')

def gene_impact_prob(muts, impacts, gene):
    """Determine the probability that a gene carrying a series of variants is neutral"""
    probs = []
    for mut_id in muts:
        imp = impacts.loc[(impacts['mut_id'] == mut_id) & (impacts['gene'] == gene)]
        if imp.type.item() == "nonsynonymous":
            probs.append(snp_neutral(mut_id, imp))

        elif imp.type.item() == "nonsense":
            probs.append(0.01)
            break

        elif imp.type.item() == "frameshift":
            probs.append(0.01)
            break

    return 1 - np.prod(probs)


def snp_neutral(mut_id, impact):
    """Determine the probability that a mutation is functionally neutral"""
    return 1

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('genotypes', metavar='G', help="Genotype table")

    parser.add_argument('mutations', metavar='M', help="Mutation table")

    parser.add_argument('--strains', '-s', default='',
                        help="List of strains (per line or comma separated)")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
