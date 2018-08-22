# Script to run intensive calculations on ko_growth data
library(tidyverse)
library(magrittr)
library(gplots)
library(ggpubr)

#### Functions ####
## Turn tibble into a matrix with rownames from a column
tbl_to_matrix <- function(x, row){
  mat <- as.matrix(select(x, -one_of(row)))
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
  mat <- as.matrix(select(x, -one_of(var)))
  rownames(mat) <- pull(x, !!var)
  cors <- cor(t(mat), method = cor_meth, use = cor_use)
  if (upper_tri){
    cors[lower.tri(cors, diag = TRUE)] <- NA
  }
  var1 <- paste0(var,'1')
  var2 <- paste0(var,'2')
  cors %<>% as_tibble(rownames=var1) %>%
    gather(key = !!var2, value = 'cor', -one_of(var1)) %>%
    drop_na()
  return(cors)
}


#### Import and Process Data ####
ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  mutate(condition = gsub('  ', ' ', condition)) %>% # Some conditions have double spaces in names
  mutate(name = if_else(is.na(name), gene, name))

strain_combs <- combn(c('S288C', 'UWOP', 'Y55', 'YPS'), 2)
strain_ko_scores <- sapply(unique(ko_growth$strain), split_strains, var='score', tbl=ko_growth, simplify = FALSE)
strain_ko_qvalues <- sapply(unique(ko_growth$strain), split_strains, var='qvalue', tbl=ko_growth, simplify = FALSE)

sig_genes <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  mutate(name = if_else(is.na(name), gene, name)) %>%
  filter(qvalue < 0.01) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
  select(-position) %>%
  unite(comb, score, qvalue) %>%
  spread(key = strain, value = comb) %>%
  mutate(num_strains = (!is.na(S288C)) + (!is.na(UWOP)) + (!is.na(Y55)) + (!is.na(YPS))) %>%
  filter(num_strains > 0) %>%
  pull(name) %>%
  unique()


#### Analysis ####
strain_gene_cors <- bind_rows(lapply(sapply(unique(ko_growth$strain), split_strains, var='score', tbl=ko_growth,
                                            col='condition', row='name', simplify = FALSE),
                                     get_cor, var='name'), .id = 'strain') %>%
  rename(gene1 = name1, gene2 = name2) %>%
  arrange(strain, gene1, gene2) %>%
  spread(strain, cor) %>%
  unite(col = 'pair', gene1, gene2, remove = FALSE)

p_strain_gene_gene_cors <- mapply(function(x, y){
  ggplot(strain_gene_cors, aes_string(x=x, y=y)) +
    geom_point(size=0.3) +
    geom_abline(intercept = 0, slope = 1, colour='firebrick2') +
    geom_abline(intercept = -0.2, slope = 1, colour='firebrick2', linetype=2) +
    geom_abline(intercept = 0.2, slope = 1, colour='firebrick2', linetype=2)
},
strain_combs[1,], strain_combs[2,], SIMPLIFY = FALSE)
names(p_strain_gene_gene_cors) <- paste(strain_combs[1,], strain_combs[2,], sep='_')

p_strain_gene_gene_cors_arr <- ggarrange(plotlist = p_strain_gene_gene_cors, ncol = 3, nrow = 2, common.legend = TRUE)
ggsave('strain_gene_gene_cors_all.jpg', p_strain_gene_gene_cors_arr, width = 7, height = 5)

p_strain_gene_gene_cors_dens <- mapply(function(x, y){
  ggplot(strain_gene_cors, aes_string(x=x, y=y)) +
    geom_bin2d(bins=40) +
    geom_abline(intercept = 0, slope = 1, colour='firebrick2') +
    geom_abline(intercept = -0.2, slope = 1, colour='firebrick2', linetype=2) +
    geom_abline(intercept = 0.2, slope = 1, colour='firebrick2', linetype=2)
},
strain_combs[1,], strain_combs[2,], SIMPLIFY = FALSE)
names(p_strain_gene_gene_cors_dens) <- paste(strain_combs[1,], strain_combs[2,], sep='_')

p_strain_gene_gene_cors_dens_arr <- ggarrange(plotlist = p_strain_gene_gene_cors_dens, ncol = 3, nrow = 2, common.legend = TRUE)
ggsave('strain_gene_gene_cors_dens_all.jpg', p_strain_gene_gene_cors_dens_arr, width = 7, height = 5)

p_strain_gene_gene_cors_contour <- mapply(function(x, y){
  ggplot(strain_gene_cors, aes_string(x=x, y=y)) +
    geom_density_2d() +
    geom_abline(intercept = 0, slope = 1, colour='firebrick2') +
    geom_abline(intercept = -0.2, slope = 1, colour='firebrick2', linetype=2) +
    geom_abline(intercept = 0.2, slope = 1, colour='firebrick2', linetype=2)
},
strain_combs[1,], strain_combs[2,], SIMPLIFY = FALSE)
names(p_strain_gene_gene_cors_contour) <- paste(strain_combs[1,], strain_combs[2,], sep='_')

p_strain_gene_gene_cors_contour_arr <- ggarrange(plotlist = p_strain_gene_gene_cors_contour, ncol = 3, nrow = 2, common.legend = TRUE)
ggsave('strain_gene_gene_cors_contour_all.jpg', p_strain_gene_gene_cors_contour_arr, width = 7, height = 5)


strain_gene_cors_mat <- filter(ko_growth, strain == 'S288C') %>%
  select(name, condition, score) %>%
  spread(key = name, value = score) %>%
  tbl_to_matrix(x = ., row = 'condition') %>%
  cor(x = ., use='pairwise')

jpeg('s288c_gene_cor_heatmap.jpg', width = 10000, height = 10000)
cols <- colorRampPalette(c('red', 'white','blue'))
heatmap.2(strain_gene_cors_mat, symm = TRUE, revC = TRUE, col = cols(200), trace = 'none')
dev.off()

ind <- rownames(strain_gene_cors_mat) %in% sig_genes
strain_gene_cors_mat_sig <- strain_gene_cors_mat[ind, ind]

jpeg('s288c_sig_gene_cor_heatmap.jpg', width = 10000, height = 10000)
cols <- colorRampPalette(c('red', 'white','blue'))
heatmap.2(strain_gene_cors_mat_sig, symm = TRUE, revC = TRUE, col = cols(200), trace = 'none')
dev.off()
