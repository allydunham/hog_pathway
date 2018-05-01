# Script to check P(Aff) scores have expected characteristics
setwd('~/Projects/hog/')
library(gplots)
library(moments)
library(tidyverse)
library(magrittr)
library(ggpubr)

#### Import Data ####
strains <- readRDS('data/Rdata/strain_meta.rds')
filtered_strains <- filter(strains, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)

genes <- readRDS('data/Rdata/gene_meta_all.rds')
hog_genes <- readRDS('data/Rdata/gene_meta_hog.rds')

complexes <- readRDS('data/Rdata/complex_members.rds')

essential <- readRDS('data/Rdata/essential_genes.rds')
essential_hash <- structure(essential$essential, names=essential$locus)

growth <- readRDS('data/Rdata/growth_liti.rds') %>%
  filter(strain %in% filtered_strains)

probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% filtered_strains, strain %in% growth$strain)

prob_mat <- readRDS('data/Rdata/paff_all_genes_mat.rds') %>%
  filter(strain %in% filtered_strains, strain %in% growth$strain)

counts <- readRDS('data/Rdata/hog_gene_mut_counts.rds') %>%
  filter(strain %in% growth$strain)

#### Analysis ####
gene_summary <- tibble(id = unique(probs$gene)) %>%
  mutate(exp_rate = apply(prob_mat, 2, function(x){fitdistr(x, densfun = 'exponential')$estimate})) %>%
  mutate(mean = colMeans(prob_mat)) %>%
  mutate(median = apply(prob_mat, 2, median)) %>%
  mutate(mean_non_zero = apply(prob_mat, 2, function(x){mean(x[!x == 0])})) %>%
  mutate(sum_zero = apply(prob_mat, 2, function(x){sum(x==0)})) %>%
  mutate(sum_non_zero = apply(prob_mat, 2, function(x){sum(!x==0)})) %>%
  left_join(genes, by = 'id') %>%
  rename(chrom=`#chrom`) %>%
  mutate(length=abs(stop - start)) %>%
  mutate(essential = essential_hash[id]) %>%
  mutate(complex = sapply(.$id, function(x){x %in% complexes$ORF}))

# Length relates to P(Aff)
p_length_mean <- ggplot(gene_summary, aes(x=length, y=mean, colour=sum_zero)) +
  geom_point() +
  xlab('Gene Length') +
  ylab('Mean P(Aff)') + 
  guides(colour=guide_colourbar(title = "Strains with P(Aff) = 0"))

ggsave('figures/paff_checks/length_vs_mean_paff.pdf', width = 12, height = 10, plot = p_length_mean)

p_length_median <- ggplot(gene_summary, aes(x=length, y=median, colour=sum_zero)) +
  geom_point() +
  xlab('Gene Length') +
  ylab('Median P(Aff)') + 
  guides(colour=guide_colourbar(title = "Strains with P(Aff) = 0"))

ggsave('figures/paff_checks/length_vs_median_paff.pdf', width = 12, height = 10, plot = p_length_median)

p_length_sum_non_zero <- ggplot(gene_summary, aes(x=length, y=sum_non_zero)) +
  geom_point() +
  xlab('Gene Length') +
  ylab('Strains where P(Aff) != 0')

ggsave('figures/paff_checks/length_vs_sum_non_zero_paff.pdf', width = 12, height = 10, plot = p_length_sum_non_zero)


probs %<>% mutate(length=structure(gene_summary$length, names=gene_summary$id)[gene]) %>%
  mutate(len_bin = cut(log10(.$length), 10))

p_length_density <- ggplot(probs, aes(x=length, y=p_aff)) +
  geom_bin2d() + scale_fill_gradient(trans="log10", guide = guide_colourbar(title = "Count")) + 
  xlab('Gene Length') + 
  ylab('P(Aff)')

ggsave('figures/paff_checks/length_vs_paff_density.pdf', width = 12, height = 10, plot = p_length_density)


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

# Essential genes
# Re-add the lower triangle of correlation matrix
gene_gene_cor[lower.tri(gene_gene_cor)] <- t(gene_gene_cor)[lower.tri(gene_gene_cor)]

# Calculate summary statistics
gene_summary %<>% mutate(cor_mean = sapply(.$id, function(x){mean(gene_gene_cor[x,], na.rm = TRUE)})) %>%
  mutate(cor_skew = sapply(.$id, function(x){skewness(gene_gene_cor[x,], na.rm = TRUE)}))

p_essential_mean <- ggplot(gene_summary, aes(x=essential, y=mean)) + geom_boxplot() + xlab('') + ylab('Mean P(Aff)')
p_essential_mean_non_zero <- ggplot(gene_summary, aes(x=essential, y=mean_non_zero)) + geom_boxplot() + xlab('') + ylab('Mean P(Aff) != 0')
p_essential_sum_zero <- ggplot(gene_summary, aes(x=essential, y=sum_zero)) + geom_boxplot() + xlab('') + ylab('Count of P(Aff) = 0')
p_essential_exp_rate <- ggplot(gene_summary, aes(x=essential, y=exp_rate)) + geom_boxplot() + ylim(0,10) + xlab('') + ylab('P(Aff) Exponential Fit Rate')
p_essential_cor_mean <- ggplot(gene_summary, aes(x=essential, y=cor_mean)) + geom_boxplot() + xlab('') + ylab('Mean P(Aff) Correlation')

p_essential <- ggarrange(p_essential_mean, p_essential_mean_non_zero, p_essential_sum_zero, p_essential_exp_rate, p_essential_cor_mean)

ggsave('figures/paff_checks/essential_genes.pdf', p_essential, width = 14, height = 10)

counts_melt <- gather(counts, key = 'gene', value = 'count', -strain) %>%
  mutate(ypdnacl1m=structure(growth$ypdnacl1m, names=growth$strain)[strain]) %>%
  mutate(ypdnacl15m=structure(growth$ypdnacl15m, names=growth$strain)[strain]) %>%
  mutate(ypdkcl2m=structure(growth$ypdkcl2m, names=growth$strain)[strain]) %>%
  mutate(ypdchx1=structure(growth$ypdchx1, names=growth$strain)[strain]) %>%
  mutate(essential=structure(essential$essential, names=essential$locus)[gene])

p_essential_mut_counts_bar <- ggplot(counts_melt, aes(x=count, y=..prop.., fill=essential)) + 
  geom_bar(position = "dodge") +
  scale_y_continuous(labels = scales::percent)

p_essential_mut_counts_box <- ggplot(counts_melt, aes(y=count, x=essential)) + 
  geom_boxplot()

p_essential_mut_counts <- ggarrange(p_essential_mut_counts_box, p_essential_mut_counts_bar)
ggsave('figures/paff_checks/essential_genes_counts.pdf', plot = p_essential_mut_counts, width = 14, height = 10)

# Complexes
p_complex_paf <- ggplot(gene_summary, aes(x=complex, y=sum_zero)) +
  geom_boxplot()

p_complex_cors_box <- ggplot(cor_genes_melt, aes(x=complex, y=cor)) +
  geom_boxplot(notch = TRUE, varwidth = TRUE) +
  xlab('Proteins in a complex?') + 
  ylab('Correlation Coefficient')

p_complex_cors_dens <- ggplot(cor_genes_melt, aes(x=cor, colour=complex)) +
  geom_density() +
  xlab('Correlation Coefficient')

ggsave('figures/paff_checks/complex_correlations_box.jpg', p_complex_cors_box, width = 12, height = 10)
ggsave('figures/paff_checks/complex_correlations_density.pdf', p_complex_cors_dens, width = 12, height = 10)


# P(Aff) and Genetic distance
genotypes <- read_tsv('data/all-genes-no-missing.genotype', col_names = TRUE)
distance_to_ref <- colSums(select(genotypes, -mut_id))

probs %<>% mutate(distance_to_ref = distance_to_ref[strain])

p_dist_ref <- ggplot(probs, aes(x=distance_to_ref, y=p_aff)) +
  geom_bin2d() +
  scale_fill_gradient(trans="log10", guide = guide_colourbar(title = "Count")) +
  xlab("Differences from SC288") + 
  ylab("P(Aff)")

ggsave('figures/paff_checks/dist_vs_p_aff_2dbin.pdf', p_dist_ref, width = 12, height = 10)

# Per strain mean P(Aff)
strain_summary <- tibble(strain = rownames(prob_mat)) %>%
  mutate(mean_paff = rowMeans(prob_mat)) %>%
  mutate(distance_to_ref = distance_to_ref[strain]) %>%
  mutate(num_high = rowSums(prob_mat > 0.9)[strain]) %>%
  mutate(num_low = rowSums(prob_mat < 0.1)[strain])

p_dist_ref_strain_mean <- ggplot(strain_summary, aes(x=distance_to_ref, y=mean_paff)) +
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean P(Aff)")

ggsave('figures/paff_checks/strain_dist_vs_paff.pdf', p_dist_ref_strain_mean, width = 12, height = 10)

p_dist_ref_strain_counts <- ggplot(strain_summary, aes(x=distance_to_ref)) +
  geom_point(aes(y=num_high, colour="P(Aff) > 0.9")) +
  geom_point(aes(y=num_low, colour="P(Aff) < 0.1")) +
  xlab("Differences from SC288") + 
  ylab("Number of Genes")

ggsave('figures/paff_checks/strain_dist_vs_paff_count_high.pdf', p_dist_ref_strain_counts, width = 12, height = 10)


# Doesn't appear to impact P(Aff) distribution
probs %<>% mutate(bin=cut(distance_to_ref, breaks = 10))

p_dist_ref_density <- ggplot(probs, aes(x=p_aff, fill=bin)) + 
  geom_density(alpha=0.5)

