#!/usr/bin/env python3
"""
Process per mutation genotypes to group them by gene and assign per gene disruption scores
If no methods are given only SIFT scores are used.

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
    impacts = pd.read_table(args.mutations, sep='\t', header=0, true_values=['True'],
                            false_values=['False'], low_memory=False)

    # Remove all synonymous entries which are not of interest
    impacts = impacts[impacts['type'] != 'synonymous']

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
    if args.conf:
        impacts['high_conf'] = ((impacts['prop_aa'] < args.conf) &
                                (impacts['type'].isin(['frameshift', 'nonsense'])))

        gene_eval_func = high_conf_gene_ko
    elif args.total:
        gene_eval_func = sum_muts
    elif args.worst:
        gene_eval_func = worst_mut
    else:
        gene_eval_func = gene_impact_prob

    if args.sift or args.blosum or args.foldx:
        snp_neutral = generate_snp_neutral(blosum=args.blosum,
                                           foldx=args.foldx,
                                           sift=args.sift)
    else:
        snp_neutral = generate_snp_neutral(blosum=False,
                                           foldx=False,
                                           sift=True)

    print("strain", *genes, sep='\t')
    for strain, genes in strain_muts.items():
        probs = {}
        for gene, mut_ids in genes.items():
            probs[gene] = gene_eval_func(muts=mut_ids,
                                         impacts=impacts,
                                         gene=gene,
                                         snp_neutral=snp_neutral)

        print(strain, *[probs[i] for i in genes], sep='\t')

def sum_muts(**kwargs):
    """Count the number of mutations in a gene"""
    return len(kwargs["muts"])

def high_conf_gene_ko(**kwargs):
    """Determine if a gene is very likely knocked out"""
    return int(kwargs["impacts"].loc[(kwargs["impacts"]['mut_id'].isin(kwargs["muts"])) &
                                     (kwargs["impacts"]['gene'] == kwargs["gene"])].high_conf.any())

def worst_mut(**kwargs):
    """Determine the probability that a gene carrying a series of variants is not
    neutral based on the highest impact variant"""
    muts = kwargs['muts']
    impacts = kwargs['impacts']
    gene = kwargs['gene']
    snp_neutral = kwargs['snp_neutral']

    # Filter impacts
    impacts = impacts[(impacts['mut_id'].isin(muts)) & (impacts['gene'] == gene)]

    # Fetch probs
    probs = get_neutral_probabilities(muts, impacts, snp_neutral)

    # Return P(Aff) = 1 - P(Neut; worst mutant)
    return 1 - min(probs)

def gene_impact_prob(**kwargs):
    """Determine the probability that a gene carrying a series of variants is not neutral"""
    muts = kwargs['muts']
    impacts = kwargs['impacts']
    gene = kwargs['gene']
    snp_neutral = kwargs['snp_neutral']

    # Filter impacts
    impacts = impacts[(impacts['mut_id'].isin(muts)) & (impacts['gene'] == gene)]

    # Fetch probs
    probs = get_neutral_probabilities(muts, impacts, snp_neutral)

    # Return P(Aff) = 1 - Prod(P(Neut))
    return 1 - np.prod(probs)

def get_neutral_probabilities(muts, impacts, snp_neutral):
    """Generate an ordered list of P(Neutral) values from a list of mut_ids in a gene"""
    # Sort mutations
    if not muts:
        return [1]

    muts = sorted(muts, key=lambda x: impacts[(impacts['mut_id'] == x)]['pos_aa'].item())

    # Calculate probabilities
    probs = []
    for mut_id in muts:
        imp = impacts.loc[(impacts['mut_id'] == mut_id)]
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
    return probs

def generate_snp_neutral(blosum=True, foldx=True, sift=True):
    """Generate a function to determine P(Neutral) for a SNP using a subset of
       the three methods"""
    funcs = []
    if blosum:
        funcs.append(lambda x: 0.66660 + 0.08293 * x.blosum62.item())

    if foldx:
        funcs.append(lambda x: 1/(1 + np.exp(0.21786182 * x.foldx_ddG.item() + 0.07351653)))

    if sift:
        funcs.append(lambda x: 1/(1 + np.exp(-1.312424 *
                                             np.log(x.sift_score.item() + 1.598027e-05) -
                                             4.103955)))

    def snp_neutral(impact):
        """Determine the probability that a mutation is functionally neutral"""
        return min(1, *[f(impact) for f in funcs])

    return snp_neutral

# def snp_neutral(impact):
#     """Determine the probability that a mutation is functionally neutral"""
#     sift = 1/(1 + np.exp(-1.312424 * np.log(impact.sift_score.item() + 1.598027e-05) - 4.103955))
#     foldx = 1/(1 + np.exp(0.21786182 * impact.foldx_ddG.item() + 0.07351653))
#     #blosum = 0.66660 + 0.08293 * impact.blosum62.item()
#     return min(1, sift, foldx) #, blosum) Better way to do blosum?

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

            if mut_id in impacts['mut_id'].values:
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

    parser.add_argument('--conf', '-c', default=0, type=float,
                        help="Only consider high confidence variants (frameshift and\
                              early stop occuring in the first X portion of the protein)")

    parser.add_argument('--total', '-t', action="store_true",
                        help="Return the count of variants in each gene")

    parser.add_argument('--worst', '-w', action="store_true",
                        help="Return the neutral probability of the most impactful variant\
                              in a gene")

    parser.add_argument('--sift', '-i', action="store_true",
                        help="Use SIFT scores to calculate P(Neutral)")

    parser.add_argument('--blosum', '-b', action="store_true",
                        help="Use BLOSUM62 scores to calculate P(Neutral)")

    parser.add_argument('--foldx', '-f', action="store_true",
                        help="Use FoldX ddG scores to calculate P(Neutral)")

    return parser.parse_args()

if __name__ == "__main__":
    ARGS = parse_args()
    main(ARGS)
