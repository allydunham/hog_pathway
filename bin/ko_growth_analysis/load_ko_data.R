# Initialisation script to load the data for KO growht analyis

# Load Packages
library(rlang)
library(tidyverse)
library(magrittr)
library(broom)
library(GSA)

## Gene Meta
gene_meta <- readRDS('data/Rdata/gene_meta_all.rds') %>%
  mutate(name = ifelse(is.na(name), id, name)) %>%
  rename(gene = id)

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

# Version of ko_growth with non sig s-scores replaced with 0 to reduce heatmap noise
ko_growth_sig <- mutate(ko_growth, score = if_else(qvalue < 0.01, score, 0))

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

set_genes <- unlist(sets) %>% unique()

set_meta <- bind_rows(sapply(sets, function(x){data_frame(set_size=length(x))}, simplify = FALSE), .id = 'gene_set')

## Load Reactome Pathways
pathways <- read_tsv('meta/reactomePathwaysScerevisiae.tsv', col_names = c('gene', 'pathID', 'http', 'path', 'IEA', 'species')) %>%
  select(gene, pathID, path) %>%
  left_join(., select(gene_meta, gene, name), by='gene') %>%
  mutate(name = ifelse(is.na(name), gene, name)) %>%
  group_by(path) %>%
  do(genes=c(.$name))
pathways <- structure(pathways$genes, names=pathways$path)
pathway_meta <- bind_rows(sapply(pathways, function(x){data_frame(set_size=length(x))}, simplify = FALSE), .id = 'gene_set')

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
growth_diff_cons <- c('Paraquat (48H)', 'Caffeine 20mM (48H)', 'Caffeine 15mM (48H)', 'Acetic acid (48H)', 'Maltose 2% (48H)', 'Maltose 2% (72H)')

## Switching probabilities (probability difference in growth between strains for a gene ko is by chance)
switch_probs <- read_tsv('data/gene_switch_probabilities.tsv') %>%
  left_join(., select(gene_meta, gene, name))

## Genotypes and mut impacts
frq <- readRDS('data/Rdata/allele_freqs.rds')
geno <- readRDS('data/Rdata/genotypes_all_genes.rds')
imp <- readRDS('data/Rdata/all_muts_impacts.rds') %>%
  left_join(., frq, by='mut_id') %>%
  left_join(., gene_meta, by = 'gene')

## RNA seq
rna_seq <- sapply(dir('data/raw/ko_BY4742/', pattern = '[0-9].csv'), function(x){read_csv(paste0('data/raw/ko_BY4742/', x), col_names = TRUE) %>%
    rename(gene = X1) %>%
    left_join(., select(gene_meta, gene, name))}, simplify = FALSE) %>%
  set_names(gsub('.csv', '', names(.))) %>%
  bind_rows(., .id = 'strain')

## Some padj found to be NA despite presence of test statistic - not unusually distributed
p <- ggplot(rna_seq, aes(x=strain, y=baseMean, col=is.na(padj))) + geom_boxplot()

rna_seq_caff <- sapply(dir('data/raw/ko_BY4742/', pattern = '.caff.csv'), function(x){read_csv(paste0('data/raw/ko_BY4742/', x), col_names = TRUE) %>%
    rename(gene = X1) %>%
    left_join(., select(gene_meta, gene, name))}, simplify = FALSE) %>%
  set_names(gsub('.csv', '', names(.))) %>%
  bind_rows(., .id = 'strain')
