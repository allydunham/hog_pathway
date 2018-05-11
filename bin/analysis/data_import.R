# Script to import and process data table TSVs into RDS files for easy analysis
# No filtering is done to strains with growth data - must be done on import
setwd('~/Projects/hog/')
library(tidyverse)
library(magrittr)

#### Meta Info ####
gene_meta_hog <- read_table2("meta/hog-gene-loci", col_names = TRUE) %>% 
  rename(chrom=`#chrom`)
saveRDS(gene_meta_hog, file = 'data/Rdata/gene_meta_hog.rds')

gene_meta_all <- read_tsv('meta/sacc_gene_loci', col_names = TRUE) %>% 
  rename(chrom=`#chrom`)
saveRDS(gene_meta_all, file = 'data/Rdata/gene_meta_all.rds')

strain_meta <- read_tsv('meta/strain_information.tsv', col_names = TRUE, comment = "") %>%
  mutate(`Isolate name`=str_to_upper(`Isolate name`))
saveRDS(strain_meta, file = 'data/Rdata/strain_meta.rds')

strain_to_sys <- structure(strain_meta$`Standardized name`, names=strain_meta$`Isolate name`)

complexes <- read_tsv('data/complexes.tsv', col_names = TRUE)
saveRDS(complexes, file = 'data/Rdata/complex_members.rds')

essential <- read_tsv('data/raw/ogee_gene_ko_lethal_350.tsv', col_names = TRUE, skip = 5) %>%
  rename(tax_D = `#taxID`)
saveRDS(essential, file = 'data/Rdata/essential_genes.rds')

#### Genotypes ####
genotypes_hog <- read_tsv("data/hog-gene-variants.all-genotypes", col_names = TRUE) %>%
  gather(strain, val, -mut_id) %>%
  spread(mut_id, val) %>%
  arrange(strain)
saveRDS(genotypes_hog, file = 'data/Rdata/genotypes_hog_genes.rds')

genotypes_all <- read_tsv('data/all-genes-no-missing.genotype', col_names = TRUE)
saveRDS(genotypes_all, file = 'data/Rdata/genotypes_all_genes.rds')

#### Variant Impacts ####
## Use filtered version to avoid poor annotations
impacts <- read_tsv('data/all-genes-no-missing.impact.filtered', col_names = TRUE, na = "NA",
                    col_types = cols(foldx_int_ddG = col_double(), foldx_int_ddG_sd = col_double()))
saveRDS(impacts, 'data/Rdata/all_muts_impacts.rds')

#### KO Probabilities ####
prob_aff_hog <- read_tsv("data/hog-gene-variants.probs.blosum", col_names = TRUE) %>%
  select(strain, gene_meta_hog$id) %>%
  set_names(c('strain', gene_meta_hog$name)) %>%
  gather(key = 'gene', value = 'p_aff', -strain)
saveRDS(prob_aff_hog, file = 'data/Rdata/paff_hog_genes.rds')

prob_aff_all <- read_tsv('data/all-genes-no-missing.koprob', col_names = TRUE) 
saveRDS(prob_aff_all, file = 'data/Rdata/paff_all_genes_mat.rds')

prob_aff_all %<>% gather(key = 'gene', value = 'p_aff', -strain)
saveRDS(prob_aff_all, file = 'data/Rdata/paff_all_genes.rds')

worst_probs_hog <- read_tsv('data/hog-gene-variants.worst-probs', col_names = TRUE)
saveRDS(worst_probs_hog, 'data/Rdata/worst_probs_hog.rds')

worst_probs_all <- read_tsv('data/all-genes-no-missing.worst-prob', col_names = TRUE)
saveRDS(worst_probs_all, 'data/Rdata/worst_probs_all_mat.rds')

worst_probs_all %<>% gather(key='gene', value = 'worst_p_aff', -strain)
saveRDS(worst_probs_all, 'data/Rdata/worst_probs_all.rds')

#### Growth Data ####
growth_bede <- read_tsv("data/raw/bede_2017_parsed.tsv", col_names = TRUE) %>%
  mutate(strain=str_to_upper(strain)) %>%
  filter(!is.na(strain_to_sys[strain])) %>%
  mutate(strain=strain_to_sys[strain])
saveRDS(growth_bede, file = 'data/Rdata/growth_bede.rds')

growth_liti <- read_tsv(file = 'data/raw/phenoMatrix_35ConditionsNormalizedByYPD.tab', col_names = TRUE) %>% 
  rename(strain=X1) %>%
  set_names(str_to_lower(names(.)))
saveRDS(growth_liti, file = 'data/Rdata/growth_liti.rds')

growth_bede_gene_dels <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(qvalue < 0.05) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  unite(comb, score, qvalue) %>%
  spread(key = strain, value = comb) %>%
  mutate(num_strains = (!is.na(S288C)) + (!is.na(UWOP)) + (!is.na(Y55)) + (!is.na(YPS))) %>%
  separate(S288C, c('S288C-score', 'S288C-qvalue'), sep = '_') %>%
  separate(UWOP, c('UWOP-score', 'UWOP-qvalue'), sep = '_') %>%
  separate(Y55, c('Y55-score', 'Y55-qvalue'), sep = '_') %>%
  separate(YPS, c('YPS-score', 'YPS-qvalue'), sep = '_') %>%
  mutate_at(vars(contains('-')), funs(as.numeric))
saveRDS(growth_bede_gene_dels, 'data/Rdata/gene_ko_growth.rds')

#### Genetic Distance ####
load('data/Rdata/genetic_distance_old.Rdata')
saveRDS(genetic_distance, 'data/Rdata/genetic_distance_matrix.rds')

genetic_distance[lower.tri(genetic_distance)] <- NA
genetic_distance_melt <- as_tibble(genetic_distance, rownames = 'strain') %>%
  gather(key = 'strain2', value = 'distance', -strain) %>%
  drop_na(distance)
saveRDS(genetic_distance_melt, 'data/Rdata/genetic_distance.rds')

saveRDS(growth_distance_NaCl6, 'data/Rdata/growth_distance_bede_NaCl6.rds')
saveRDS(growth_distance_NaCl4, 'data/Rdata/growth_distance_bede_NaCl4.rds')

#### Pathway Data ####
hog_path <- read_tsv(file = 'data/hog-gene-variants.path.blosum', col_names = TRUE)
  mutate(hog_active, hog_active = as.factor(hog_active))
saveRDS(hog_path, 'data/Rdata/hog_path_probs.rds')
  
hog_counts <- read_tsv('data/hog-gene-variants.mut-counts', col_names = TRUE)
saveRDS(hog_counts, 'data/Rdata/hog_gene_mut_counts.rds')

all_counts <- read_tsv('data/all-genes-no-missing.mut-counts', col_names = TRUE)
saveRDS(all_counts, 'data/Rdata/all_gene_mut_counts.rds')

#### Correlation Data ####
load('data/correlation_data.Rdata')
gene_gene_cor[lower.tri(gene_gene_cor)] <- t(gene_gene_cor)[lower.tri(gene_gene_cor)]
diag(gene_gene_cor) <- 1
saveRDS(gene_gene_cor, 'data/Rdata/all_gene_correlations_matrix.rds')
saveRDS(cor_genes_melt, 'data/Rdata/all_gene_correlations.rds')
saveRDS(cor_growth, 'data/Rdata/gene_growth_correlations_matrix.rds')
saveRDS(cor_growth_melt, 'data/Rdata/gene_growth_correlations.rds')

#### Frequencies ####
allele_freqs <- read_tsv('data/all-genes-no-missing.mut-freqs', col_names = TRUE)
saveRDS(allele_freqs, 'data/Rdata/allele_freqs.rds')

#### Normalised P(Aff) Data ####
probs_norm <- read_tsv('data/all-genes-no-missing.koprob.norm', col_names = TRUE)
saveRDS(probs_norm, 'data/Rdata/norm_paff_mat.rds')
probs_norm %<>% gather(key = 'gene', value = 'norm_p_aff', -strain)
saveRDS(probs_norm, 'data/Rdata/norm_paff.rds')

counts_norm <- read_tsv('data/all-genes-no-missing.mut-counts.norm', col_names = TRUE)
saveRDS(counts_norm, 'data/Rdata/norm_counts_mat.rds')
counts_norm %<>% gather(key = 'gene', value = 'norm_count', -strain)
saveRDS(counts_norm, 'data/Rdata/norm_counts.rds')

worst_norm <- read_tsv('data/all-genes-no-missing.worst-prob.norm', col_names = TRUE)
saveRDS(worst_norm, 'data/Rdata/norm_worst_mat.rds')
worst_norm %<>% gather(key = 'gene', value = 'norm_worst_p_aff', -strain)
saveRDS(worst_norm, 'data/Rdata/norm_worst.rds')

