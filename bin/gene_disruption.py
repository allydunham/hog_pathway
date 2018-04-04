#!/usr/bin/env python3
"""
Process per mutation genotypes to group them by gene and assign per gene disruption scores

ToDo
- Currently nothing takes account of heterozygotes
- blosum integration
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
    # Import impact table
    impacts = pd.read_table(args.mutations)

    # Extract list of affected genes
    genes = impacts.gene.dropna().unique()

    # Import genotype data
    strain_muts = import_genotypes(args.genotypes, genes, impacts, args.strains)

    # Print set of files with mutations in each gene - dir: args.print/strain/gene_file
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

    # Determine per gene mutation probabilities and print table
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
    #blosum = 0.66660 + 0.08293 * impact.blosum62.item()
    return min(1, sift, foldx) #, blosum) Better way to do blosum?

def import_genotypes(genotype_file, genes, impacts, strains=''):
    """Import a table of genotypes into a nested dictionary
       With optional strain filtering"""
    # Import genotypes and split into genes
    with fileinput.input(genotype_file) as geno_file:
        header = next(geno_file).strip().split('\t')

        # Read in strains to filter
        if strains:
            try:
                strains = [line.strip() for line in fileinput.input(strains)]
            except FileNotFoundError:
                strains = strains.strip().split(',')
        else:
            strains = header[1:]

        # Get column indeces for desired strains
        indeces = [i for i in range(len(header)) if header[i] in strains]

        # Set up storage dict - strain_muts[strain][gene][muts]
        strain_muts = {s:{g:[] for g in genes} for s in strains}

        # Extract mutations per gene per strain
        for line in geno_file:
            line = line.strip().split('\t')
            mut_id = line[0]

            # Identify impacted genes
            affected_genes = impacts[impacts['mut_id'] == mut_id].gene.dropna().values

            # Add mutations to appropriate strain/gene
            for i in indeces:
                if not line[i] == '0':
                    for gene in affected_genes:
                        strain_muts[header[i]][gene].append(mut_id)

    return strain_muts

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
