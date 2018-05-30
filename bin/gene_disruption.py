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
import sys
import numpy as np
import pandas as pd

def main(args):
    """Main script"""
    # Import impact table
    impacts = pd.read_table(args.mutations, sep='\t', header=0, true_values=['True'],
                            false_values=['False'], low_memory=False)

    # Remove all synonymous entries which are not of interest
    impacts = impacts[impacts['type'] != 'synonymous']

    # Determine SNP evaluation function and apply it to variants
    if args.sift or args.blosum or args.foldx:
        snp_neutral = generate_snp_neutral(blosum=args.blosum,
                                           foldx=args.foldx,
                                           sift=args.sift)
    else:
        snp_neutral = generate_snp_neutral(blosum=False,
                                           foldx=True,
                                           sift=True)

    impacts['p_neut'] = impacts.apply(snp_neutral, axis=1)

    # Sort by order variants occur in genes
    impacts = impacts.sort_values(by='prop_aa', axis=0)

    # Extract list of affected genes and generate a default P(Aff) series
    genes = impacts.gene.dropna().unique()

    # Import genotype data
    strain_muts = pd.read_table(args.genotypes, sep='\t', header=0,
                                low_memory=False, index_col=0)

    # Determine function to evaluate genes
    if args.conf:
        # Only count high confidence variants and give the total in a gene
        impacts['high_conf'] = ((impacts['prop_aa'] < args.conf) &
                                (impacts['type'].isin(['frameshift', 'nonsense'])))
        gene_eval_func = lambda x: x['high_conf'].sum()
        process_nonsense = True
    elif args.total:
        # Give the total number of variants
        gene_eval_func = lambda x: x.shape[0]
        process_nonsense = False
    elif args.worst:
        gene_eval_func = lambda x: 1 - x['p_neut'].min()
        process_nonsense = True
    else:
        gene_eval_func = lambda x: 1 - x['p_neut'].prod()
        process_nonsense = True

    strain_probs = strain_muts.apply(func=eval_strain, axis=0, reduce=False,
                                     impacts=impacts, genes=genes,
                                     eval_func=gene_eval_func,
                                     process_nonsense=process_nonsense,
                                     zygosity=args.zygosity).transpose()
    strain_probs.index.name = 'strain'
    strain_probs.to_csv(sys.stdout, sep='\t')

def eval_strain(muts, impacts, genes, eval_func, process_nonsense=True, zygosity='both'):
    """
    Evaluate a strain genotype to give a list of gene P(Aff) values
    """
    # Set default P(Aff) as 0
    default_probs = pd.Series([0] * len(genes), index=genes)

    # Extract variants in strain
    if zygosity == 'hom':
        mut_ids = muts[muts > 1].index
    elif zygosity == 'het':
        mut_ids = muts[muts == 1].index
    else:
        mut_ids = muts[muts > 0].index

    # Determine P(Aff) for each affected gene
    gene_probs = impacts[(impacts['mut_id'].isin(mut_ids))
                        ].groupby('gene').apply(gene_impact_prob,
                                                eval_func=eval_func,
                                                process_nonsense=process_nonsense)
    default_probs[gene_probs.index] = gene_probs

    return default_probs

def gene_impact_prob(gene_imp, eval_func, process_nonsense=True):
    """Determine the P(Aff) for a gene carrying a table of variants"""
    # Filter to first nonsense/framshift
    if process_nonsense and ('nonsense' in gene_imp['type'] or
                             'frameshift' in gene_imp['type']):
        first_index = gene_imp[gene_imp['type'].isin(('nonsense', 'frameshift'))].index[0]
        gene_imp = gene_imp.iloc[:gene_imp.index.get_loc(first_index) + 1,]

    # Combine probabilities
    return eval_func(gene_imp)

def generate_snp_neutral(blosum=True, foldx=True, sift=True):
    """Generate a function to determine P(Neutral) for a SNP using a subset of
       the three methods"""
    funcs = []
    if blosum:
        funcs.append(lambda x: 0.66660 + 0.08293 * x['blosum62'])

    if foldx:
        funcs.append(lambda x: 1/(1 + np.exp(0.21786182 * x['foldx_ddG'] + 0.07351653)))

    if sift:
        funcs.append(lambda x: 1/(1 + np.exp(-1.312424 *
                                             np.log(x['sift_score'] + 1.598027e-05) -
                                             4.103955)))

    def snp_neutral(impact):
        """Determine the probability that a mutation is functionally neutral"""
        if impact['type'] == "nonsynonymous":
            return min(1, *[f(impact) for f in funcs])

        elif impact['type'] in ('nonsense', 'frameshift'):
            if impact['prop_aa'] < 0.95:
                return 0.01

            return 0.99

        return 1

    return snp_neutral

def parse_args():
    """Process input arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('genotypes', metavar='G', help="Genotype table")

    parser.add_argument('mutations', metavar='M', help="Mutation table")

    parser.add_argument('--zygosity', '-z', default='both', type=str,
                        help="Minimum number of variant alleles required to be variant\
                              (1=het, 2=hom for diploids)")

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
