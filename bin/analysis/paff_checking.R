# Script to check P(Aff) scores have expected characteristics
setwd('~/Projects/hog/')
library(gplots)
library(MASS)
library(moments)
library(tidyverse)
library(magrittr)
library(ggpubr)

#### Import Data ####
strains <- readRDS('data/Rdata/strain_meta.rds')
filtered_strains <- filter(strains, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)
# Filter strains that are unusually distant from the reference genome
filtered_strains <-  setdiff(filtered_strains , c("AMH", "BAG", "BAH", "BAL", "CEG", "CEI"))

genes <- readRDS('data/Rdata/gene_meta_all.rds')
hog_genes <- readRDS('data/Rdata/gene_meta_hog.rds')

complexes <- readRDS('data/Rdata/complex_members.rds')

essential <- readRDS('data/Rdata/essential_genes.rds')
essential_hash <- structure(essential$essential, names=essential$locus)

impacts <- readRDS('data/Rdata/all_muts_impacts.rds')

growth <- readRDS('data/Rdata/growth_liti.rds') %>%
  filter(strain %in% filtered_strains)

genotypes <- readRDS('data/Rdata/genotypes_all_genes.rds')
distance_to_ref <- structure(strains$`Total number of SNPs`, names=strains$`Standardized name`)

probs <- readRDS('data/Rdata/paff_all_genes.rds') %>%
  filter(strain %in% filtered_strains, strain %in% growth$strain) %>%
  mutate(distance_to_ref = distance_to_ref[strain]) %>% 
  mutate(length=structure(genes$stop - genes$start, names=genes$id)[gene]) %>%
  mutate(len_bin = cut(log10(.$length), 10))

prob_mat <- readRDS('data/Rdata/paff_all_genes_mat.rds') %>%
  filter(strain %in% filtered_strains, strain %in% growth$strain) %>%
  select(-strain) %>%
  as.matrix() %>%
  set_rownames(growth$strain)  

worst_probs <- readRDS('data/Rdata/worst_probs_all.rds') %>%
  filter(strain %in% intersect(filtered_strains, growth$strain))

counts <- readRDS('data/Rdata/all_gene_mut_counts.rds') %>%
  filter(strain %in% growth$strain)

allele_freqs <- readRDS('data/Rdata/allele_freqs.rds')

gene_gene_cor <- readRDS('data/Rdata/all_gene_correlations_matrix.rds')

#### Analysis ####
gene_summary <- tibble(id = unique(probs$gene)) %>%
  mutate(exp_rate = apply(prob_mat, 2, function(x){fitdistr(x, densfun = 'exponential')$estimate})) %>%
  mutate(mean = colMeans(prob_mat)) %>%
  mutate(median = apply(prob_mat, 2, median)) %>%
  mutate(mean_non_zero = apply(prob_mat, 2, function(x){mean(x[!x == 0])})) %>%
  mutate(sum_zero = apply(prob_mat, 2, function(x){sum(x==0)})) %>%
  mutate(sum_non_zero = apply(prob_mat, 2, function(x){sum(!x==0)})) %>%
  left_join(genes, by = 'id') %>%
  mutate(length=abs(stop - start)) %>%
  mutate(essential = essential_hash[id]) %>%
  mutate(complex = sapply(.$id, function(x){x %in% complexes$ORF}))

# Length relates to P(Aff)
fit_length <- lm(mean ~ length, data = gene_summary)
probs %<>% mutate(length_corrected = p_aff - (
  fit_length$coefficients[1] + fit_length$coefficients[2] * length))
gene_summary %<>% mutate(mean_length_corrected = mean - (fit_length$coefficients[1] + fit_length$coefficients[2] * length))

p_length_mean <- ggplot(gene_summary, aes(x=length, y=mean, colour=sum_zero)) +
  geom_point() +
  xlab('Gene Length') +
  ylab('Mean P(Aff)') + 
  guides(colour=guide_colourbar(title = "Strains with P(Aff) = 0")) + 
  stat_function(fun = function(x){fit_length$coefficients[1] + fit_length$coefficients[2] * x})

ggsave('figures/paff_checks/length_vs_mean_paff.pdf', width = 12, height = 10, plot = p_length_mean)

p_length_mean_corrected <- ggplot(gene_summary, aes(x=length, y=mean_length_corrected, colour=sum_zero)) +
  geom_point() +
  xlab('Gene Length') +
  ylab('Mean Length Corrected P(Aff)') + 
  guides(colour=guide_colourbar(title = "Strains with P(Aff) = 0")) + 
  geom_smooth(method = 'lm')

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

## Essential genes
# Calculate summary statistics
gene_summary %<>% mutate(cor_mean = colMeans(gene_gene_cor, na.rm = TRUE)[id]) %>%
  mutate(cor_skew = apply(gene_gene_cor, 2, skewness,  na.rm = TRUE)[id])

p_essential_mean <- ggplot(gene_summary, aes(x=essential, y=mean)) + geom_boxplot() + xlab('') + ylab('Mean P(Aff)')
p_essential_mean_non_zero <- ggplot(gene_summary, aes(x=essential, y=mean_non_zero)) + geom_boxplot() + xlab('') + ylab('Mean P(Aff) != 0')
p_essential_sum_zero <- ggplot(gene_summary, aes(x=essential, y=sum_zero)) + geom_boxplot() + xlab('') + ylab('Count of P(Aff) = 0')
p_essential_exp_rate <- ggplot(gene_summary, aes(x=essential, y=exp_rate)) + geom_boxplot() + ylim(0,10) + xlab('') + ylab('P(Aff) Exponential Fit Rate')
p_essential_cor_mean <- ggplot(gene_summary, aes(x=essential, y=cor_mean)) + geom_boxplot() + xlab('') + ylab('Mean P(Aff) Correlation')

p_essential <- ggarrange(p_essential_mean, p_essential_mean_non_zero, p_essential_sum_zero, p_essential_exp_rate, p_essential_cor_mean)

ggsave('figures/paff_checks/essential_genes.pdf', p_essential, width = 14, height = 10)

# Gene variant counts
counts_melt <- gather(counts, key = 'gene', value = 'count', -strain)

probs %<>% mutate(variant_count = counts_melt$count) %>%
  mutate(ypdnacl1m=unname(structure(growth$ypdnacl1m, names=growth$strain)[strain])) %>%
  mutate(ypdnacl15m=unname(structure(growth$ypdnacl15m, names=growth$strain)[strain])) %>%
  mutate(ypdkcl2m=unname(structure(growth$ypdkcl2m, names=growth$strain)[strain])) %>%
  mutate(ypdchx1=unname(structure(growth$ypdchx1, names=growth$strain)[strain])) %>%
  mutate(essential=NA)

# Assign essentiallity (structure method is too slow)
essential_genes <- filter(essential, essential == 'E') %>% pull(locus)
non_essential_genes <- filter(essential, essential == 'NE') %>% pull(locus)
probs[probs$gene %in% essential_genes, 'essential'] <- 'E'
probs[probs$gene %in% non_essential_genes, 'essential'] <- 'NE'

p_essential_mut_counts_bar <- ggplot(probs, aes(x=count, y=..prop.., fill=essential)) + 
  geom_bar(position = "dodge") +
  scale_y_continuous(labels = scales::percent) %>%
  xlab('Numbr of Variants') %>%
  ylab('Proportion of Genes')

p_essential_mut_counts_box <- ggplot(probs, aes(y=count, x=essential)) + 
  geom_boxplot() +
  scale_y_log10()

p_essential_mut_counts <- ggarrange(p_essential_mut_counts_box, p_essential_mut_counts_bar)
ggsave('figures/paff_checks/essential_genes_counts.pdf', plot = p_essential_mut_counts, width = 14, height = 10)

# Different but not a meaningful amount
wilcox.test(filter(counts_melt, essential == 'E') %>% pull(count), filter(counts_melt, essential == 'NE') %>% pull(count))

# Most impactful variant
probs %<>% mutate(worst_p_aff = worst_probs$worst_p_aff)

p_worst_paff_essential_box <- ggplot(filter(probs, !is.na(essential), variant_count > 0), aes(x = essential, y = worst_p_aff)) +
  geom_boxplot()

p_worst_paff_essential_density <- ggplot(filter(probs, !is.na(essential), variant_count > 0), aes(colour = essential, x = worst_p_aff)) +
  geom_density()

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
p_dist_ref <- ggplot(probs, aes(x=distance_to_ref, y=p_aff)) +
  geom_bin2d() +
  scale_fill_gradient(trans="log10", guide = guide_colourbar(title = "Count")) +
  xlab("Differences from SC288") + 
  ylab("P(Aff)")

ggsave('figures/paff_checks/dist_vs_p_aff_2dbin.pdf', p_dist_ref, width = 12, height = 10)

# Per strain mean P(Aff)
strain_summary <- tibble(strain = rownames(prob_mat)) %>%
  mutate(mean_paff = rowMeans(prob_mat)) %>%
  mutate(median_paff = apply(prob_mat, 1, median)) %>%
  mutate(distance_to_ref = distance_to_ref[strain]) %>%
  mutate(num_high = rowSums(prob_mat > 0.9)[strain]) %>%
  mutate(num_low = rowSums(prob_mat < 0.1)[strain])

fit_strain <- lm(mean_paff ~ distance_to_ref, data = strain_summary)
strain_summary %<>% mutate(mean_corrected = mean_paff - (fit_strain$coefficients[1] + fit_strain$coefficients[2] * distance_to_ref))

p_dist_ref_strain_mean <- ggplot(strain_summary, aes(x=distance_to_ref, y=mean_paff)) +
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean P(Aff)") + 
  geom_smooth(method = 'lm')

ggsave('figures/paff_checks/strain_dist_vs_mean_paff.pdf', p_dist_ref_strain_mean, width = 12, height = 10)

p_dist_ref_strain_mean_corrected <- ggplot(strain_summary, aes(x=distance_to_ref, y=mean_corrected)) +
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Mean Phylogeny Corrected P(Aff)") + 
  geom_smooth(method = 'lm')

p_dist_ref_strain_counts <- ggplot(strain_summary, aes(x=distance_to_ref)) +
  geom_point(aes(y=num_high, colour="P(Aff) > 0.9")) +
  geom_point(aes(y=num_low, colour="P(Aff) < 0.1")) +
  xlab("Differences from SC288") + 
  ylab("Number of Genes")

ggsave('figures/paff_checks/strain_dist_vs_paff_count_high.pdf', p_dist_ref_strain_counts, width = 12, height = 10)

p_dist_ref_strain_median <- ggplot(strain_summary, aes(x=distance_to_ref, y=median_paff)) +
  geom_point() +
  xlab("Differences from SC288") + 
  ylab("Median P(Aff)")

ggsave('figures/paff_checks/strain_dist_vs_median_paff.pdf', p_dist_ref_strain_median, width = 12, height = 10)

# Doesn't appear to impact P(Aff) distribution
probs %<>% mutate(bin=cut(distance_to_ref, breaks = 10))

p_dist_ref_density <- ggplot(probs, aes(x=p_aff, fill=bin)) + 
  geom_density(alpha=0.5)

#### Analyse Individual Variants ####
impacts %<>% mutate(freq=structure(allele_freqs$freq, names=allele_freqs$mut_id)[mut_id]) %>%
  mutate(p_neut_sift= 1/(1 + exp(-1.312424 * log(sift_score + 1.598027e-05) - 4.103955))) %>%
  mutate(p_neut_foldx= 1/(1 + exp(0.21786182 * foldx_ddG + 0.07351653))) %>%
  mutate(p_neut_blosum = 0.66660 + 0.08293 * blosum62) %>%
  filter(!type == 'synonymous') %>%
  mutate(essential = essential_hash[gene])

p_freq_neut <- ggplot(impacts, aes(x=freq)) +
  geom_point(aes(y=p_neut_sift, colour="SIFT")) +
  geom_point(aes(y=p_neut_foldx, colour="FoldX")) +
  geom_point(aes(y=p_neut_blosum, colour="BLOSUM62")) + 
  xlab('Allele Frequency') + 
  ylab('P(Neutral)') + 
  guides(colour=guide_legend(title = 'Method'))
ggsave('figures/paff_checks/freq_vs_p_neut.pdf', p_freq_neut, width = 12, height = 10)  

p_freq_essential <- ggplot(impacts, aes(x=essential, y=freq)) + 
  geom_boxplot() + 
  scale_y_log10()
ggsave('figures/paff_checks/freq_vs_essential.pdf', p_freq_essential, width = 12, height = 10)  

p_freq_essential_sift <- ggplot(filter(impacts, sift_score < 0.05, !is.na(essential)), aes(x=essential, y=freq)) + 
  geom_boxplot() + 
  scale_y_log10()

p_freq_essential_foldX <- ggplot(filter(impacts, abs(foldx_ddG) > 2, !is.na(essential)), aes(x=essential, y=freq)) + 
  geom_boxplot() + 
  scale_y_log10()

# But does not appear significant
t.test(freq ~ essential, data = filter(impacts, sift_score < 0.05, !is.na(essential)))
t.test(freq ~ essential, data = filter(impacts, abs(foldx_ddG) > 2, !is.na(essential)))

p_essential_foldx <- ggplot(filter(impacts, !is.na(essential)), aes(colour=essential, x=sift_score)) +
  geom_density()

impact_summary <- group_by(impacts, gene) %>%
  summarise(essential = first(essential),
            count_sift = sum(sift_score < 0.05, na.rm = TRUE),
            count_foldx = sum(foldx_ddG > 2, na.rm = TRUE)) %>%
  mutate(length = structure(gene_summary$length, names=gene_summary$id)[gene]) %>%
  mutate(count_sift_per_base = count_sift / length) %>%
  mutate(count_foldx_per_base = count_foldx / length)

# Observe same result as Omars Mutfunc paper
p_essential_sift_per_base <- ggplot(filter(impact_summary, !is.na(essential)), aes(x=essential, y=count_sift_per_base)) + 
  geom_boxplot(notch = TRUE, varwidth = TRUE) +
  xlab('Gene Essential?') +
  ylab('Number of variants with SIFT < 0.05 per base')

t.test(count_sift_per_base ~ essential, data = impact_summary)

ggsave('figures/paff_checks/gene_count_low_sift_per_base_essential_box.pdf', p_essential_sift_per_base, width = 12, height = 10)

# but less so here?
p_essential_foldx_per_base <- ggplot(impact_summary, aes(x=essential, y=count_foldx_per_base)) + 
  geom_boxplot(notch = TRUE, varwidth = TRUE)

t.test(count_foldx_per_base ~ essential, data = impact_summary)

# Integrate frequency of variants into analysis?



