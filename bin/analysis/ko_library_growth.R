# Script performing analying ko library growth
setwd('~/Projects/hog/')
library(rlang)
library(tidyverse)
library(magrittr)
library(ggpubr)
library(gplots)
library(ggdendro)
library(dendextend)
library(plotly)

#### Import and Process Data ####
genes <- readRDS('data/Rdata/gene_meta_all.rds')
essential <- readRDS('data/Rdata/essential_genes.rds')
essential_genes <- filter(essential, essential == 'E') %>% pull(locus)
non_essential_genes <- filter(essential, essential == 'NE') %>% pull(locus)

ko_growth <- read_tsv('data/raw/ko_scores.txt', col_names = TRUE) %>%
  filter(!gene == 'WT') %>%
  filter(!duplicated(.[,c('strain', 'condition', 'gene')])) %>%
#  filter(qvalue < 0.5) %>% # Possibly worthwhile fitlering to only important kos?
  select(-position) %>%
  mutate(condition = gsub('  ', ' ', condition)) # Some conditions have double spaces in names

strain_combs <- combn(c('S288C', 'UWOP', 'Y55', 'YPS'), 2)

# Does not apply sorting because ko_scores tsv comes sorted
split_strains <- function(str, var, tbl){
  # Function to filter ko_growth to a given strain with a matrix of condition against 'var' (score or qvalue)
  ko <- filter(tbl, strain == str) %>%
    select(condition, gene, !!var) %>%
    spread(key = gene, value = !!var)
  return(ko)
}

strain_ko_scores <- lapply(unique(ko_growth$strain), split_strains, var='score', tbl=ko_growth)
names(strain_ko_scores) <- unique(ko_growth$strain)

strain_ko_qvalues <- lapply(unique(ko_growth$strain), split_strains, var='qvalue', tbl=ko_growth)
names(strain_ko_qvalues) <- unique(ko_growth$strain)

#### Analysis ####
## Dendograms
get_dend <- function(x){
  mat <- as.matrix(select(x, -condition))
  rownames(mat) <- x$condition
  return(as.dendrogram(hclust(dist(mat))))
}

strain_con_dends <- lapply(strain_ko_scores, get_dend)

p_dend_list <- lapply(strain_con_dends, ggdendrogram, rotate=TRUE)
p_dends <- ggarrange(plotlist = p_dend_list, ncol = 2, nrow = 2, labels = names(strain_con_dends))
ggsave('figures/ko_growth/strain_condition_dendograms.pdf', p_dends, width = 10, height = 10)

# dendogram comparisons
# must filter Maltose 2%  (72H) since no data is available for Y55
strain_ko_scores_filtered <- lapply(unique(ko_growth$strain), split_strains, var='score', tbl=filter(ko_growth, !condition == 'Maltose 2% (72H)'))
names(strain_ko_scores_filtered) <- unique(ko_growth$strain)
strain_con_dends_filtered <- lapply(strain_ko_scores_filtered, get_dend)

pdf('figures/ko_growth/strain_condition_comparison_dendograms.pdf', width = 15, height = 20)
layout(matrix(1:(6*3), nrow = 3, byrow = TRUE), widths = rep(c(2,0.75,2), 2))
for (c in 1:dim(strain_combs)[2]){
  strains <- strain_combs[,c]
  tanglegram(untangle(strain_con_dends_filtered[[strains[1]]], strain_con_dends_filtered[[strains[2]]]),
             margin_inner = 15, common_subtrees_color_branches = TRUE,
             main_left = strains[1], main_right = strains[2], lwd=2, edge.lwd=2, just_one = FALSE)
}
layout(matrix(1, nrow = 1))
dev.off()


## Correlations
get_cor <- function(x, cor_meth = 'pearson', cor_use = 'pairwise', upper_tri=TRUE){
  mat <- as.matrix(select(x, -condition))
  rownames(mat) <- x$condition
  cors <- cor(t(mat), method = cor_meth, use = cor_use)
  if (upper_tri){
    cors[lower.tri(cors, diag = TRUE)] <- NA
  }
  cors %<>% as_tibble(rownames='condition1') %>%
    gather(key = 'condition2', value = 'cor', -condition1) %>%
    drop_na()
  return(cors)
}

strain_con_cors <- bind_rows(lapply(strain_ko_scores, get_cor), .id = 'strain') %>%
  arrange(strain, condition1, condition2) %>%
  spread(strain, cor) %>%
  unite(col = 'pair', condition1, condition2, remove = FALSE)

p_strain_con_cor_cors <- mapply(function(x, y){
  ggplot(strain_con_cors, aes_string(x=x, y=y, text='pair')) + 
    geom_point(size=0.5) + 
    geom_abline(intercept = 0, slope = 1, colour='firebrick2') +
    geom_abline(intercept = -0.2, slope = 1, colour='firebrick2', linetype=2) +
    geom_abline(intercept = 0.2, slope = 1, colour='firebrick2', linetype=2)
  }, 
  strain_combs[1,], strain_combs[2,], SIMPLIFY = FALSE)
names(p_strain_con_cor_cors) <- paste(strain_combs[1,], strain_combs[2,], sep='_')

# removed rows correspond to pairs including Maltose 2%  (72H), since no data is available for Y55
p_strain_con_cor_cors_arr <- ggarrange(plotlist = p_strain_con_cor_cors, ncol = 3, nrow = 2, common.legend = TRUE)
ggsave('figures/ko_growth/strain_condition_condition_cors.pdf', p_strain_con_cor_cors_arr, width = 7, height = 5)

## Interactive plots
# Single
ggplotly(p_strain_con_cor_cors$UWOP_YPS, tooltip = c('pair'))

# Multiple
plotly_strain_cors <- subplot(sapply(p_strain_con_cor_cors, function(x){ggplotly(x, tooltip = c('pair'))}, simplify = FALSE),
                              nrows = 2, titleX = TRUE, titleY = TRUE, margin = c(0.05, 0.05, 0.07, 0.07))

# 3D Plots - show some joint
p_strain_cors_3d <- plot_ly(strain_con_cors, x=~S288C, y=~UWOP, z=~Y55, mode = 'markers', text=~pair) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'S288C'),
                      yaxis = list(title = 'UWOP'),
                      zaxis = list(title = 'Y55')))



## Per gene correlation
ko_growth_gene <- mutate(ko_growth, name = if_else(is.na(name),gene,name)) %>%
  select(-qvalue, -gene) %>%
  rename(gene = name) %>%
  spread(key = strain, value = score) %>%
  group_by(gene) %>%
  summarise(UWOP_S288C = cor(UWOP, S288C, use = 'na.or.complete'),
            UWOP_Y55 = cor(UWOP, Y55, use = 'na.or.complete'),
            UWOP_YPS = cor(UWOP, YPS, use = 'na.or.complete'),
            S288C_Y55 = cor(S288C, Y55, use = 'na.or.complete'),
            S288C_YPS = cor(S288C, YPS, use = 'na.or.complete'),
            YPS_Y55 = cor(YPS, Y55, use = 'na.or.complete')) %>%
  mutate(mean_cor = (UWOP_S288C + UWOP_Y55 + UWOP_YPS + S288C_Y55 + S288C_YPS + YPS_Y55)/6, 
         var_cor = (UWOP_S288C^2 + UWOP_Y55^2 + UWOP_YPS^2 + S288C_Y55^2 + S288C_YPS^2 + YPS_Y55^2)/6 - mean_cor^2)
  
p_gene_strain_cor_boxes <- ggplot(gather(ko_growth_gene, key = 'strain_pair', value = 'cor', -gene), aes(x=strain_pair, y=cor, text=gene)) +
  geom_jitter()
ggplotly(p_gene_strain_cor_boxes, tooltip = c('gene'))

p_gene_mean_var_cor <- ggplot(ko_growth_gene, aes(x=var_cor, y=mean_cor, label=gene)) +
  geom_point() +
  xlab('Variance of gene KO s-score correlation accross strains') +
  ylab('Mean gene KO s-score correlation accross strains')
ggplotly(p_gene_mean_var_cor)
