# Script to analyse whole genome P(Aff) scores and their correlations
setwd('~/Projects/hog/')
library(MASS)
library(mixtools)
library(gplots)
library(tidyverse)
library(magrittr)

#### Import ####
## Strain information
meta <- read_tsv('meta/strain_information.tsv', col_names = TRUE)


## Gene information
genes <- read_tsv('meta/sacc_gene_loci', col_names = TRUE)
sys_to_gene <- structure(genes$name, names=genes$id)
filtered_strains <- filter(meta, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)

hog_genes <- read_table2('meta/hog-gene-loci', col_names = FALSE, comment = '#') %>%
  set_names(c('chrom', 'start', 'stop', 'id', 'name', 'strand'))


# Growth data
growth <- read_tsv(file = 'data/raw/phenoMatrix_35ConditionsNormalizedByYPD.tab', col_names = TRUE) %>% 
  rename(strain=X1) %>%
  set_names(str_to_lower(names(.)))

## P(aff)
probs <- read_tsv('data/all-genes-no-missing.koprob', col_names = TRUE) %>%
  filter(strain %in% growth$strain)

# Filter genes with no variation
p_aff_sums <- colSums(select(probs, -strain))
no_prob_genes <- names(which(p_aff_sums == 0))
probs %<>% select(-one_of(no_prob_genes))

# Complex Membership
complexes <- read_tsv('data/complexes.tsv', col_names = TRUE)

#### Analysis ####
## Gene/Gene Correlation
cor_genes <- cor(select(probs, -strain))

cor_genes[lower.tri(cor_genes, diag = TRUE)] <- NA

same_complex <- function(x, y, comp=complexes){
  x_comps <- comp[comp$ORF == x,] %>% pull(Complex)
  y_comps <- comp[comp$ORF == y,] %>% pull(Complex)
  return(any(x_comps %in% y_comps))
}

cor_genes_melt <- as_tibble(cor_genes) %>%
  add_column(gene = rownames(cor_genes), .before = 1) %>%
  gather(key = 'gene2', value = 'cor', -gene) %>%
  drop_na(cor) %>%
  mutate(mag = abs(cor)) %>%
  arrange(desc(mag)) %>%
  mutate(name = sys_to_gene[gene], name2 = sys_to_gene[gene2]) %>%
  mutate(complex = map2(.$gene, .$gene2, same_complex))

# pdf('figures/all_genes_paff_heatmap.pdf', width = 50, height = 50)
# cols <- colorRampPalette(c("blue", "white","red"))(256)
# heatmap.2(cor_genes, symm = TRUE, revC = TRUE, col=cols,
#           breaks=seq(-1,1,2/256), trace = "none")
# dev.off()

## Gene/Growth Correlation
cor_growth <- cor(select(probs, -strain), select(growth, -strain))
meanGrowthCor <- rowMeans(cor_growth)

cor_growth_melt <- as_tibble(cor_growth, rownames = 'gene_id') %>%
  gather(key='condition', value = 'correlation', -gene_id) %>%
  group_by(condition) %>%
  mutate(gene_name = sys_to_gene[gene_id]) %>%
  mutate(mag = abs(correlation)) %>% 
  mutate(cor_sub_mean = correlation - meanGrowthCor) %>%
  mutate(cor_div_mean = correlation/meanGrowthCor) %>%
  arrange(desc(mag), .by_group = TRUE)

pdf('figures/all_genes_growth_heatmap.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth, col=cols, breaks=seq(-1,1,2/256), trace = "none", symkey = FALSE)
dev.off()

# normalise by subtraction
pdf('figures/all_genes_growth_heatmap_norm_sub.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth - meanGrowthCor, col=cols, breaks=seq(-1,1,2/256), trace = "none", symkey = FALSE)
dev.off()

# Normalise by division
pdf('figures/all_genes_growth_heatmap_norm_div.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth / meanGrowthCor, col=cols, trace = "none", symkey = FALSE)
dev.off()

# Hog Genes Only
pdf('figures/hog_genes_growth_heatmap.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth[rownames(cor_growth) %in% hog_genes$id,], breaks=seq(-1,1,2/256), 
          col=cols, trace = "none", symkey = FALSE, cexRow = 2.5, cexCol = 2.5, margins = c(24,15))
dev.off()


  

# Individual Gene Analysis
gene_p_aff_means <- colMeans(select(probs, -strain))

prob_melt <- gather(probs, key = 'gene', value = 'p_aff', -strain)

prob_mat <- as.matrix(select(probs, -strain))

gene_summary <- tibble(id = names(select(probs, -strain))) %>%
  mutate(exp_rate = apply(select(probs, -strain), 2, function(x){fitdistr(x, densfun = 'exponential')$estimate})) %>%
  mutate(mean = colMeans(select(probs, -strain))) %>%
  mutate(mean_non_zero = apply(prob_mat, 2, function(x){mean(x[!x == 0])})) %>%
  mutate(sum_zero = apply(prob_mat, 2, function(x){sum(x==0)})) %>%
  mutate(sum_non_zero = apply(prob_mat, 2, function(x){sum(!x==0)})) %>%
  left_join(genes, by = 'id') %>%
  rename(chrom=`#chrom`) %>%
  mutate(length=abs(stop - start))


# Length does relate to P(Aff)
p_length <- ggplot(gene_summary, aes(x=length, y=mean_non_zero)) +
  geom_point() +
  xlab('Gene Length') +
  ylab('Mean P(Aff)')
  
ggsave('figures/correlations/length_vs_mean_paff.pdf', width = 12, height = 10, plot = p_length)

# Strand is unrelated
p_strand <- ggplot(gene_summary, aes(x=strand, y=mean)) +
  geom_boxplot() +
  xlab('Srand') + 
  ylab('Mean P(Aff)')

#ks.test(filter(gene_summary, strand == '+') %>% pull(mean),
#        filter(gene_summary, strand == '-') %>% pull(mean))
# Two-sample Kolmogorov-Smirnov test
# 
# data:  filter(gene_summary, strand == "+") %>% pull(mean) and filter(gene_summary, strand == "-") %>% pull(mean)
# D = 0.021794, p-value = 0.5312
# alternative hypothesis: two-sided

# No apparant relationship with chromosome
p_chrom <- ggplot(gene_summary, aes(x=chrom, y=mean)) +
  geom_boxplot() +
  xlab('Chromosome') + 
  ylab('Mean P(Aff)') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

chrom_test <- t(sapply(unique(gene_summary$chrom), function(x){
  t <- ks.test(filter(gene_summary, chrom == x) %>% pull(mean),
              filter(gene_summary, !chrom == x) %>% pull(mean));
  return(structure(c(x, t$p.value),
                   names = c('chrom', 'p.val')))
})) %>% as_tibble() %>% mutate(p.adj = p.adjust(.$p.val, method = 'fdr'))

# Relationship between mean with and without zero entries
p_zero_non_zero <- ggplot(gene_summary, aes(x=mean, y=mean_non_zero, colour=sum_zero)) + 
  geom_point() + 
  xlab('Mean P(Aff)') + 
  ylab('Mean P(Aff) (Non-zero samples only)') + 
  coord_equal() + 
  xlim(0,1) + 
  ylim(0,1)

# Distribution of zeros
p_zero_dist <- ggplot(gene_summary, aes(x=sum_zero)) +
  geom_histogram(binwidth = 20) +
  xlab('Strains where P(Aff) = 0')

