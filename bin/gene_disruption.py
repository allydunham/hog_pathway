#!/usr/bin/env python3
"""
Process per mutation genotypes to group them by gene and assign per gene disruption scores

ToDo
- Currently nothing takes account of heterozygotes
- Option to print mutations in each gene in each strain (before/instead of collating them to probs)
- Add blosum where values are missing
- deal with start codons
"""
import argparse
import fileinput
import sys
import os
import numpy as np
import pandas as pd

def main(args):
    """Main script"""
    impacts = pd.read_table(args.mutations)

    genes = impacts.gene.dropna().unique()

    # Read in and filter strains
    if args.strains:
        try:
            strains = []
            with fileinput.input(args.strains) as strain_file:
                for line in strain_file:
                    strains.append(line.strip())

        except FileNotFoundError:
            strains = args.strains.strip().split(',')

    with fileinput.input(args.genotypes) as geno_file:
        header = next(geno_file).strip().split('\t')

        if not args.strains:
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

    if args.print:
        for strain, genes in strain_muts.items():
            for gene, muts in genes.items():
                path = ''.join((args.print, strain, '/'))
                os.makedirs(path, exist_ok=True)
                with open(''.join((path, gene, '.tsv')), 'w') as gene_file:
                    impacts[(impacts['mut_id'].isin(muts)) &
                            (impacts['gene'] == gene)].to_csv(gene_file,
                                                              sep='\t',
                                                              na_rep='NA')

        sys.exit(0)

    # Determine per gene mutation probabilities
    print("strain", *genes, sep='\t')
    for strain, genes in strain_muts.items():
        probs = {}
        for gene, muts in genes.items():
            probs[gene] = gene_impact_prob(muts, impacts, gene)

        print(strain, *[probs[i] for i in genes], sep='\t')

def gene_impact_prob(muts, impacts, gene):
    """Determine the probability that a gene carrying a series of variants is neutral"""
    # Sort mutations
    muts = sorted(muts, key=lambda x: impacts[(impacts['mut_id'] == x) &
                                              (impacts['gene'] == gene)]['pos_aa'].item())

    # Calculate probabilities
    probs = []
    for mut_id in muts:
        imp = impacts.loc[(impacts['mut_id'] == mut_id) & (impacts['gene'] == gene)]
        if imp.type.item() == "nonsynonymous":
            probs.append(snp_neutral(imp))

        elif imp.type.item() == "nonsense":
            if imp.prop_aa.item() < 0.95:
                probs.append(0.01)
            else:
                probs.append(0.99)
            break

        elif imp.type.item() == "frameshift":
            if imp.prop_aa.item() < 0.95:
                probs.append(0.01)
            else:
                probs.append(0.99)
            break
    return 1 - np.prod(probs)


def snp_neutral(impact):
    """Determine the probability that a mutation is functionally neutral"""
    sift = 1/(1 + np.exp(-1.312424 * np.log(impact.sift_score.item() + 1.598027e-05) - 4.103955))
    foldx = 1/(1 + np.exp(0.21786182 * impact.foldx_ddG.item() + 0.07351653))
    return min(1, sift, foldx)

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('genotypes', metavar='G', help="Genotype table")

    parser.add_argument('mutations', metavar='M', help="Mutation table")

    parser.add_argument('--strains', '-s', default='',
                        help="List of strains (per line or comma separated)")

    parser.add_argument('--print', '-p', default='',
                        help="Print the impact of variants in each gene per\
                              species in a specified folder")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
