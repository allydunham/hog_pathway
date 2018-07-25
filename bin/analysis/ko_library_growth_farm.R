# Script to run intensive calculations on ko_growth data
library(tidyverse)
library(magrittr)
library(gplots)

#### Functions ####
## Turn tibble into a matrix with rownames from a column
tbl_to_matrix <- function(x, row){
  mat <- as.matrix(select(x, -row))
  rownames(mat) <- pull(x, !!row)
  return(mat)
}


## Split ko_growth tibble into subtable of a strain giving a variable spread over gene or condition
# Does not apply sorting because ko_scores tsv comes sorted
split_strains <- function(str, var, tbl, col='gene', row='condition'){
  # Function to filter ko_growth to a given strain with a matrix of condition against 'var' (score or qvalue)
  ko <- filter(tbl, strain == str) %>%
    select(!!row, !!col, !!var) %>%
    spread(key = col, value = !!var)
  return(ko)
}


## Calculate correlation between pairs of genes or conditions (based on profile according to the other)
# requires a tibble with the dependant variable in rows labeled by column 'var'
get_cor <- function(x, cor_meth = 'pearson', cor_use = 'pairwise', upper_tri=TRUE, var='condition'){
  mat <- as.matrix(select(x, -!!var))
  rownames(mat) <- pull(x, !!var)
  cors <- cor(t(mat), method = cor_meth, use = cor_use)
  if (upper_tri){
    cors[lower.tri(cors, diag = TRUE)] <- NA
  }
  var1 <- paste0(var,'1')
  var2 <- paste0(var,'2')
  cors %<>% as_tibble(rownames=var1) %>%
    gather(key = !!var2, value = 'cor', -!!var1) %>%
    drop_na()
  return(cors)
}


#### Import and Process Data ####
genes <- readRDS('data/Rdata/gene_meta_all.rds')
essential <- readRDS('data/Rdata/essential_genes.rds')
essential_genes <- filter(essential, essential == 'E') %>% pull(locus)
non_essential_genes <- filter(essential, essential == 'NE') %>% pull(locus)

ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  mutate(condition = gsub('  ', ' ', condition)) %>% # Some conditions have double spaces in names
  mutate(name = if_else(is.na(name), gene, name))

strain_combs <- combn(c('S288C', 'UWOP', 'Y55', 'YPS'), 2)
strain_ko_scores <- sapply(unique(ko_growth$strain), split_strains, var='score', tbl=ko_growth, simplify = FALSE)
strain_ko_qvalues <- sapply(unique(ko_growth$strain), split_strains, var='qvalue', tbl=ko_growth, simplify = FALSE)


#### Analysis ####
strain_gene_cors <- bind_rows(lapply(sapply(unique(ko_growth$strain), split_strains, var='score', tbl=ko_growth, 
                                            col='condition', row='name', simplify = FALSE), 
                                     get_cor, var='name'), .id = 'strain') %>%
  rename(gene1 = name1, gene2 = name2) %>%
  arrange(strain, gene1, gene2) %>%
  spread(strain, cor) %>%
  unite(col = 'pair', gene1, gene2, remove = FALSE)

p_strain_gene_gene_cors <- mapply(function(x, y){
  ggplot(strain_gene_cors, aes_string(x=x, y=y, text='pair')) + 
    geom_point(size=0.5) + 
    geom_abline(intercept = 0, slope = 1, colour='firebrick2') +
    geom_abline(intercept = -0.2, slope = 1, colour='firebrick2', linetype=2) +
    geom_abline(intercept = 0.2, slope = 1, colour='firebrick2', linetype=2)
}, 
strain_combs[1,], strain_combs[2,], SIMPLIFY = FALSE)
names(p_strain_gene_gene_cors) <- paste(strain_combs[1,], strain_combs[2,], sep='_')

p_strain_gene_gene_cors_arr <- ggarrange(plotlist = p_strain_gene_gene_cors, ncol = 3, nrow = 2, common.legend = TRUE)
ggsave('strain_gene_gene_cors_all.jpg', p_strain_gene_gene_cors_arr, width = 7, height = 5)

p_strain_gene_gene_cors_dens <- mapply(function(x, y){
  ggplot(strain_gene_cors, aes_string(x=x, y=y, text='pair')) + 
    geom_bin2d() + 
    geom_abline(intercept = 0, slope = 1, colour='firebrick2') +
    geom_abline(intercept = -0.2, slope = 1, colour='firebrick2', linetype=2) +
    geom_abline(intercept = 0.2, slope = 1, colour='firebrick2', linetype=2)
}, 
strain_combs[1,], strain_combs[2,], SIMPLIFY = FALSE)
names(p_strain_gene_gene_cors_dens) <- paste(strain_combs[1,], strain_combs[2,], sep='_')

p_strain_gene_gene_cors_dens_arr <- ggarrange(plotlist = p_strain_gene_gene_cors_dens, ncol = 3, nrow = 2, common.legend = TRUE)
ggsave('strain_gene_gene_cors_dens_all.jpg', p_strain_gene_gene_cors_dens_arr, width = 7, height = 5)

strain_gene_cors_mat <- select(strain_gene_cors, gene1, gene2, S288C) %>%
  spread(key = 'gene2', value = 'S288C') %>%
  tbl_to_matrix(., 'gene1')

jpeg('figures/ko_growth/s288c_gene_cor_heatmap.jpg', width = 1000, height = 1000)
heatmap.2(strain_gene_cors_mat, symm = TRUE, revC = TRUE)
dev.off()
