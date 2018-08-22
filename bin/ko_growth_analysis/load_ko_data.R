# Initialisation script to load the data for KO growht analyis
setwd('~/Projects/hog/')
source('bin/general_functions.R')
source('bin/ko_growth_analysis/ko_analysis_functions.R')

# Load Packages
library(rlang)
library(tidyverse)
library(magrittr)
library(broom)
library(GSA)

## KO Growth data
ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  mutate(condition = gsub('  ', ' ', condition)) %>% # Some conditions have double spaces in names
  mutate(name = if_else(is.na(name), gene, name))

ko_growth_spread <- filter(ko_growth, qvalue < 0.01) %>%
  select(strain, condition, name, score) %>%
  spread(key = strain, value = score)

ko_growth_spread_all <- ko_growth %>%
  select(strain, condition, name, score) %>%
  spread(key = strain, value = score)

## Identify Tested Genes/Strains/Conditions
conditions <- unique(ko_growth$condition)
condition_combs <- combn(conditions, 2)

# Based on knowledge plus correlation of gene profiles accross different strains
equiv_cons <- structure(c("2,4-Dichlorophenoxyacetic acid", "39ºC", "39ºC", "5-FU", "6-AU", "39ºC", "39ºC", "Acetic acid",
                          "Amphotericin B", "Amphotericin B", "Anaerobic", "Cadmium chloride", "Caffeine", "Caffeine",
                          "Caspofungin", "Caspofungin", "Clozapine", "Cyclohexamide", "DMSO", "Doxorubicin", "Glucose",
                          "Glycerol", "Glycerol", "Maltose", "Maltose", "NaCl", "NaCl", "NaCl", "NaCl", "NaCl", "NaCl", "NaCl",
                          "NiSO4", "Nitrogen starvation", "Nystatin", "Paraquat", "Paraquat", "SC + hepes", "Sorbitol", 
                          "aa starvation"),
                        names = conditions)

strains <- c('S288C', 'UWOP', 'Y55', 'YPS')
strain_combs <- combn(strains, 2)

genes <- unique(ko_growth$name)

## Load Meta lists of essential genes, complexes and gene sets
essential <- readRDS('data/Rdata/essential_genes.rds')
essential_genes <- filter(essential, essential == 'E') %>% pull(locus)
non_essential_genes <- filter(essential, essential == 'NE') %>% pull(locus)

complexes <- readRDS('data/Rdata/complex_members.rds') %>%
  group_by(Complex) %>%
  do(gene=c(.$Name)) %>%
  mutate(size = length(gene),
         num_tested = sum(gene %in% genes))

# Sets split into GO types from which they were derived
# Gene sets sourced from http://www.go2msig.org/cgi-bin/prebuilt.cgi?taxid=559292
gene_sets_bp <- GSA.read.gmt('meta/cerevisiae_bp_go_gene_sets.gmt')
names(gene_sets_bp$genesets) <- gene_sets_bp$geneset.names
gene_sets_bp_filt <- gene_sets_bp$genesets[sapply(gene_sets_bp$genesets, function(x){sum(x %in% genes) > 5})]

gene_sets_mf <- GSA.read.gmt('meta/cerevisiae_mf_go_gene_sets.gmt')
names(gene_sets_mf$genesets) <- gene_sets_mf$geneset.names

gene_sets_cc <- GSA.read.gmt('meta/cerevisiae_cc_go_gene_sets.gmt')
names(gene_sets_cc$genesets) <- gene_sets_cc$geneset.names

# Combined list
gene_sets <- list(bp=gene_sets_bp$genesets, cc=gene_sets_cc$genesets, mf=gene_sets_mf$genesets)
gene_sets_filt <- lapply(gene_sets, function(sets){sets[sapply(sets, function(x){sum(x %in% genes) > 5})]})

# Flattened combined list
sets <- unlist(gene_sets_filt, recursive = FALSE)


## Extract S-Scores and Q-Values into per strain matrices
strain_ko_scores <- sapply(unique(ko_growth$strain), split_strains, var='score', tbl=ko_growth, simplify = FALSE)
strain_ko_qvalues <- sapply(unique(ko_growth$strain), split_strains, var='qvalue', tbl=ko_growth, simplify = FALSE)

## Find genes significnat in at least n conditions
sig_genes <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  mutate(name = if_else(is.na(name), gene, name)) %>%
  filter(qvalue < 0.01) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  unite(comb, score, qvalue) %>%
  spread(key = strain, value = comb) %>%
  mutate(num_strains = (!is.na(S288C)) + (!is.na(UWOP)) + (!is.na(Y55)) + (!is.na(YPS)))

sig_genes_4 <- sig_genes %>%
  filter(num_strains > 3) %>%
  pull(name) %>%
  unique()

sig_genes_3 <- sig_genes %>%
  filter(num_strains > 2) %>%
  pull(name) %>%
  unique()

sig_genes_2 <- sig_genes %>%
  filter(num_strains > 1) %>%
  pull(name) %>%
  unique()

sig_genes_1 <- sig_genes %>%
  filter(num_strains > 0) %>%
  pull(name) %>%
  unique()

## Strain/Condition KS tests for Gene sets
ks_batches <- read_tsv('data/ko_ks_tests/all_tests.tsv')

## Growth rate of each strain in each condition relative to S288C
strain_relative_growth <- read_tsv('data/relative_fitness.tsv')
