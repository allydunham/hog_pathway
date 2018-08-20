## Script to perform batch KS tests on gene sets by strain and condition
library(tidyverse)
library(magrittr)
library(broom)
library(GSA)

#### Process Arguments and set Params ####
args <- commandArgs(trailingOnly = TRUE)
print(args)
## Get list of conditions to test in batch from args
condition_inds <- args

## Number of genes in a set required to have been tested overall
minimum_genes <- 10

## Number of genes in set required in strain/condition pair to perform ks test
minimum_genes_test <- 5

#### Import Data ####
# ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
#   filter(!gene == 'WT') %>%
#   filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
#   select(-position) %>%
#   mutate(condition = gsub('  ', ' ', condition)) %>% # Some conditions have double spaces in names
#   mutate(name = if_else(is.na(name), gene, name))
# 
# genes <- unique(ko_growth$name)
# strains <- c('S288C', 'UWOP', 'Y55', 'YPS')
# conditions <- sort(unique(ko_growth$condition))
# 
# ## Gene sets
# gene_sets_bp <- GSA.read.gmt('meta/cerevisiae_bp_go_gene_sets.gmt')
# names(gene_sets_bp$genesets) <- gene_sets_bp$geneset.names
# 
# gene_sets_mf <- GSA.read.gmt('meta/cerevisiae_mf_go_gene_sets.gmt')
# names(gene_sets_mf$genesets) <- gene_sets_mf$geneset.names
# 
# gene_sets_cc <- GSA.read.gmt('meta/cerevisiae_cc_go_gene_sets.gmt')
# names(gene_sets_cc$genesets) <- gene_sets_cc$geneset.names
# 
# gene_sets <- list(bp=gene_sets_bp$genesets, cc=gene_sets_cc$genesets, mf=gene_sets_mf$genesets)
# gene_sets_filt <- lapply(gene_sets, function(sets){sets[sapply(sets, function(x){sum(x %in% genes) > minimum_genes})]})
# save.image('data/Rdata/batch_ks_tests.Rdata')

load('data/Rdata/batch_ks_tests.Rdata')
conditions <- conditions[condition_inds]
ko_growth <- filter(ko_growth, condition %in% conditions)

#### Functions ####
## Perform ks test on a gene sets vs all other genes in a given table of ko scores
set_ks_test <- function(tbl, set){
  ind <- tbl$name %in% set
  # Return NAs if not enough genes in the set have been tested
  if (sum(ind) < minimum_genes_test){
    return(data_frame(statistic=NA, p.value=NA, method='Two-sample Kolmogorov-Smirnov test', alternative='two-sided'))
  } else {
    return(tidy(ks.test(tbl[ind,]$score, tbl[!ind,]$score)))
  }
}

## Manage gene sets and labelling. Send each set in each group for ks testing sequentially and bind results
do_set_tests <- function(tbl, gene_sets){
  return(bind_rows(
          sapply(gene_sets, function(set_group){
            bind_rows(
              sapply(set_group, function(set){set_ks_test(tbl, set)}, simplify = FALSE),
              .id = 'gene_set'
            )},
            simplify = FALSE),
          .id = 'gene_set_group')
        )
}

#### Perform KS tests ####
ks_tests <- group_by(ko_growth, condition, strain) %>%
  do(do_set_tests(., gene_sets_filt))

write_tsv(ks_tests, paste0('data/ko_ks_tests/', paste0(gsub('º', '', gsub(' ','_',conditions)), collapse = '_')))
