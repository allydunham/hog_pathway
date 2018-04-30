# Script to analyse whole genome P(Aff) scores and their correlations
setwd('~/Projects/hog/')
library(MASS)
library(mixtools)
library(gplots)
library(tidyverse)
library(magrittr)
library(moments)
library(ggpubr)

# #### Import ####
# ## Strain information
# meta <- read_tsv('meta/strain_information.tsv', col_names = TRUE)
# 
# 
# ## Gene information
# genes <- read_tsv('meta/sacc_gene_loci', col_names = TRUE)
# sys_to_gene <- structure(genes$name, names=genes$id)
# filtered_strains <- filter(meta, Ploidy == 2, Aneuploidies == 'euploid') %>% pull(`Standardized name`)
# 
# hog_genes <- read_table2('meta/hog-gene-loci', col_names = FALSE, comment = '#') %>%
#   set_names(c('chrom', 'start', 'stop', 'id', 'name', 'strand'))
# 
# 
# # Growth data
# growth <- read_tsv(file = 'data/raw/phenoMatrix_35ConditionsNormalizedByYPD.tab', col_names = TRUE) %>% 
#   rename(strain=X1) %>%
#   set_names(str_to_lower(names(.)))
# 
# ## P(aff)
# probs <- read_tsv('data/all-genes-no-missing.koprob', col_names = TRUE) %>%
#   filter(strain %in% growth$strain)
# 
# # Filter genes with no variation
# p_aff_sums <- colSums(select(probs, -strain))
# no_prob_genes <- names(which(p_aff_sums == 0))
# probs %<>% select(-one_of(no_prob_genes))
# 
# # Complex Membership
# complexes <- read_tsv('data/complexes.tsv', col_names = TRUE)
# 
# ## Calculations
# # Gene/Gene
# cor_genes <- cor(select(probs, -strain))
# 
# cor_genes[lower.tri(cor_genes, diag = TRUE)] <- NA
# 
# same_complex <- function(x, y, comp=complexes){
#   x_comps <- comp[comp$ORF == x,] %>% pull(Complex)
#   y_comps <- comp[comp$ORF == y,] %>% pull(Complex)
#   return(any(x_comps %in% y_comps))
# }
# 
# cor_genes_melt <- as_tibble(cor_genes) %>%
#   add_column(gene = rownames(cor_genes), .before = 1) %>%
#   gather(key = 'gene2', value = 'cor', -gene) %>%
#   drop_na(cor) %>%
#   mutate(mag = abs(cor)) %>%
#   arrange(desc(mag)) %>%
#   mutate(name = sys_to_gene[gene], name2 = sys_to_gene[gene2]) %>%
#   mutate(complex = map2(.$gene, .$gene2, same_complex))
# 
# # Gene/Growth
# cor_growth <- cor(select(probs, -strain), select(growth, -strain))
# meanGrowthCor <- rowMeans(cor_growth)
# 
# cor_growth_melt <- as_tibble(cor_growth, rownames = 'gene_id') %>%
#   gather(key='condition', value = 'correlation', -gene_id) %>%
#   group_by(condition) %>%
#   mutate(gene_name = sys_to_gene[gene_id]) %>%
#   mutate(mag = abs(correlation)) %>% 
#   mutate(cor_sub_mean = correlation - meanGrowthCor) %>%
#   mutate(cor_div_mean = correlation/meanGrowthCor) %>%
#   arrange(desc(mag), .by_group = TRUE)

load('data/correlation_data.Rdata')

hog_genes <- read_table2('meta/hog-gene-loci', col_names = FALSE, comment = '#') %>%
  set_names(c('chrom', 'start', 'stop', 'id', 'name', 'strand'))

counts <- read_tsv('data/hog-gene-variants.mut-counts', col_names = TRUE) %>%
  filter(strain %in% growth$strain)

## Import gene essentiallity
essential <- read_tsv('data/raw/ogee_gene_ko_lethal_350.tsv', col_names = TRUE, skip = 5) %>%
  rename(tax_D = `#taxID`)

essential_hash <- structure(essential$essential, names=essential$locus)

## Probs melt
prob_mat <- as.matrix(select(probs, -strain)) %>% set_rownames(probs$strain)
probs %<>% gather(key = 'gene', value = 'p_aff', -strain)

#### Analysis ####
## Gene/Gene Correlation
pdf('figures/heatmaps/all_genes_paff_heatmap.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_genes, symm = TRUE, revC = TRUE, col=cols,
          breaks=seq(-1,1,2/256), trace = "none")
dev.off()

## Gene/Growth Correlation
pdf('figures/heatmaps/all_genes_growth_heatmap.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth, col=cols, breaks=seq(-1,1,2/256), trace = "none", symkey = FALSE)
dev.off()

# normalise by subtraction
pdf('figures/heatmaps/all_genes_growth_heatmap_norm_sub.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth - meanGrowthCor, col=cols, breaks=seq(-1,1,2/256), trace = "none", symkey = FALSE)
dev.off()

# Normalise by division
pdf('figures/heatmaps/all_genes_growth_heatmap_norm_div.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth / meanGrowthCor, col=cols, trace = "none", symkey = FALSE)
dev.off()

# Hog Genes Only
cor_hog_counts <- cor(select(counts, -strain))
rownames(cor_hog_counts) <- structure(hog_genes$name, names=hog_genes$id)[rownames(cor_hog_counts)]

pdf('figures/heatmaps/hog_genes_mut_count_heatmap.pdf', width = 20, height = 20)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_hog_counts, symm = TRUE, revC = TRUE, col=cols,
          breaks=seq(-1,1,2/256), trace = "none")
dev.off()

# growth
cor_growth_hog <- cor_growth[rownames(cor_growth) %in% hog_genes$id,]
rownames(cor_growth_hog) <- structure(hog_genes$name, names=hog_genes$id)[rownames(cor_growth_hog)]
  
pdf('figures/heatmaps/hog_genes_growth_heatmap.pdf', width = 50, height = 50)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth_hog, breaks=seq(-1,1,2/256), 
          col=cols, trace = "none", symkey = FALSE, cexRow = 2.5, cexCol = 2.5, margins = c(24,15))
dev.off()

# Only select conditions
cor_growth_hog_osmotic <- cor_growth_hog[,c("ypdkcl2m","ypdnacl1m", "ypdnacl15m", "ypd14",
                                            "ypd40", "ypdcuso410mm", "ypdsodiummetaarsenite")]

pdf('figures/heatmaps/hog_genes_growth_osmotic_cons_heatmap.pdf', width = 15, height = 25)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_growth_hog_osmotic, breaks=seq(-1,1,2/256), 
          col=cols, trace = "none", symkey = FALSE, cexRow = 2.5, cexCol = 2.5, margins = c(30,15),
          labCol = c("KCl 2M", "NaCl 1M", "NaCl 1.5M", expression("14"~degree~"C"), expression("40"~degree~"C"), 
                                      expression("CuSO"[4]~10~"mM"), "Na Metaarsenite 2.5mM"))
dev.off()

# Against counts
cor_hog_growth_counts <- cor(select(counts, -strain), select(growth, -strain))
rownames(cor_hog_growth_counts) <- structure(hog_genes$name, names=hog_genes$id)[rownames(cor_hog_growth_counts)]

pdf('figures/heatmaps/hog_genes_mut_count_growth_heatmap.pdf', width = 20, height = 20)
cols <- colorRampPalette(c("blue", "white","red"))(256)
heatmap.2(cor_hog_growth_counts, col=cols, symkey=FALSE, cexRow = 2.5, cexCol = 2.5, margins = c(30,15),
          breaks=seq(-1,1,2/256), trace = "none")
dev.off()

# Individual Gene Analysis

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

# Length does relate to P(Aff)
p_length_mean <- ggplot(gene_summary, aes(x=length, y=mean, colour=sum_zero)) +
  geom_point() +
  xlab('Gene Length') +
  ylab('Mean P(Aff)') + 
  guides(colour=guide_colourbar(title = "Strains with P(Aff) = 0"))
  
ggsave('figures/correlations/length_vs_mean_paff.pdf', width = 12, height = 10, plot = p_length_mean)

p_length_median <- ggplot(gene_summary, aes(x=length, y=median, colour=sum_zero)) +
  geom_point() +
  xlab('Gene Length') +
  ylab('Median P(Aff)') + 
  guides(colour=guide_colourbar(title = "Strains with P(Aff) = 0"))

ggsave('figures/correlations/length_vs_median_paff.pdf', width = 12, height = 10, plot = p_length_median)

p_length_sum_non_zero <- ggplot(gene_summary, aes(x=length, y=sum_non_zero)) +
  geom_point() +
  xlab('Gene Length') +
  ylab('Strains where P(Aff) != 0')

ggsave('figures/correlations/length_vs_sum_non_zero_paff.pdf', width = 12, height = 10, plot = p_length_sum_non_zero)


probs %<>% mutate(length=structure(gene_summary$length, names=gene_summary$id)[gene]) %>%
  mutate(len_bin = cut(log10(.$length), 10))

p_length_density <- ggplot(probs, aes(x=length, y=p_aff)) +
  geom_bin2d() + scale_fill_gradient(trans="log10", guide = guide_colourbar(title = "Count")) + 
  xlab('Gene Length') + 
  ylab('P(Aff)')

ggsave('figures/correlations/length_vs_paff_density.pdf', width = 12, height = 10, plot = p_length_density)


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

ggsave('figures/correlations/essential_genes.pdf', p_essential, width = 14, height = 10)

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
ggsave('figures/correlations/essential_genes_counts.pdf', plot = p_essential_mut_counts, width = 14, height = 10)

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

ggsave('figures/correlations/complex_correlations_box.jpg', p_complex_cors_box, width = 12, height = 10)
ggsave('figures/correlations/complex_correlations_density.pdf', p_complex_cors_dens, width = 12, height = 10)


# P(Aff) and Genetic distance
genotypes <- read_tsv('data/all-genes-no-missing.genotype', col_names = TRUE)
distance_to_ref <- colSums(select(genotypes, -mut_id))

probs %<>% mutate(distance_to_ref = distance_to_ref[strain])

p_dist_ref <- ggplot(probs, aes(x=distance_to_ref, y=p_aff)) +
  geom_bin2d() +
  scale_fill_gradient(trans="log10", guide = guide_colourbar(title = "Count")) +
  xlab("Differences from SC288") + 
  ylab("P(Aff)")
  
ggsave('figures/correlations/dist_vs_p_aff_2dbin.pdf', p_dist_ref, width = 12, height = 10)

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

ggsave('figures/correlations/strain_dist_vs_paff.pdf', p_dist_ref_strain_mean, width = 12, height = 10)

p_dist_ref_strain_counts <- ggplot(strain_summary, aes(x=distance_to_ref)) +
  geom_point(aes(y=num_high, colour="P(Aff) > 0.9")) +
  geom_point(aes(y=num_low, colour="P(Aff) < 0.1")) +
  xlab("Differences from SC288") + 
  ylab("Number of Genes")

ggsave('figures/correlations/strain_dist_vs_paff_count_high.pdf', p_dist_ref_strain_counts, width = 12, height = 10)


# Doesn't appear to impact P(Aff) distribution
probs %<>% mutate(bin=cut(distance_to_ref, breaks = 10))

p_dist_ref_density <- ggplot(probs, aes(x=p_aff, fill=bin)) + 
  geom_density(alpha=0.5)


