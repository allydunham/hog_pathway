setwd("~/Projects/hog/")
library(Biostrings)

genotypes <- t(read.table("data/hog-gene-variants.all-genotypes", header=TRUE, row.names = 1))

growth <- read.table("data/raw/bede_2017_parsed.tsv",header=TRUE,sep='\t')

meta <- read.table('meta/strain_information.tsv', sep = '\t', header=TRUE, fill=TRUE, quote = "", comment.char = '')
strain_to_sys <- meta$Standardized.name
names(strain_to_sys) <- meta$Isolate.name
sys_to_strain <- names(strain_to_sys)
names(sys_to_strain) <- strain_to_sys

growth <- growth[!is.na(strain_to_sys[rownames(growth)]),]
rownames(growth) <- strain_to_sys[rownames(growth)]

genotypes <- genotypes[rownames(growth),]

genetic_distance <- as.matrix(dist(genotypes, method = "man"))
growth_distance_NaCl4 <- as.matrix(dist(growth[,"sodium.chloride.0.4mM"]))
growth_distance_NaCl6 <- as.matrix(dist(growth[,"sodium.chloride.0.6mM"]))
growth_distance_Sorb <- as.matrix(dist(growth[,"Sorbitol.1mM"]))

plot(genetic_distance[upper.tri(genetic_distance)], growth_distance_NaCl4[upper.tri(growth_distance_NaCl4)],
     pch=20, xlab="Genetic Distance (Manhatten)", ylab="Difference in S-Score (NaCl 0.4mM)",
     ylim=c(0,7))

plot(genetic_distance[upper.tri(genetic_distance)], growth_distance_NaCl6[upper.tri(growth_distance_NaCl6)],
     pch=20, xlab="Genetic Distance (Manhatten)", ylab="Difference in S-Score (NaCl 0.6mM)",
     ylim=c(0,7))

plot(genetic_distance[upper.tri(genetic_distance)], growth_distance_Sorb[upper.tri(growth_distance_Sorb)],
     pch=20, xlab="Genetic Distance (Manhatten)", ylab="Difference in S-Score (Sorbitol 1mM)",
     ylim=c(0,7))

geno_pca <- prcomp(genotypes)


#### Whole genome genetic distance ####
load('data/genetic_distance.Rdata')
library(gplots)
cols <- colorRampPalette(c("white","red"))(256)

rownames(genetic_distance) <- sys_to_strain[rownames(genetic_distance)]

pdf('figures/genetic_dist_heatmap_small.pdf', width = 20, height = 20)
heatmap.2(genetic_distance, symm = TRUE, revC = TRUE, col=cols,
          breaks=seq(0,max(genetic_distance),max(genetic_distance)/256), trace = "none")
dev.off()

# Genetic distance vs growth across all conds
library(ggplot2)
library(reshape2)

# Extract distance matrix across all conditions
growth_distance <- lapply(colnames(growth), function(x){as.matrix(dist(growth[,x]))})
names(growth_distance) <- colnames(growth)

for (i in 1:length(growth_distance)){
  dimnames(growth_distance[[i]]) <- list(rownames(growth), rownames(growth))
}

growth_melt <- melt(as.matrix(growth), variable.name = 'Condition', value.name = 'SScore')

# Set up data frame for all pairs over all conditions (excluding matching strains)
gen_growth_dist <- genetic_distance_growth
gen_growth_dist[lower.tri(gen_growth_dist, diag = TRUE)] <- NA
gen_growth_dist <- melt(gen_growth_dist, varnames = c('Strain1', 'Strain2'))
gen_growth_dist <- subset(gen_growth_dist, !is.na(value))

colnames(gen_growth_dist) <- c('Strain1', 'Strain2', 'GeneticDistance')

for (i in names(growth_distance)){
  t <- growth_distance[[i]]
  t[lower.tri(t, diag = TRUE)] <- NA
  t <- subset(melt(t), !is.na(value))
  gen_growth_dist[, i] <- t$value
}

gen_growth_dist_melt <- melt(gen_growth_dist, id.vars = c('Strain1', 'Strain2', 'GeneticDistance'),
                        variable.name = 'Condition', value.name = 'SScoreDiff')

p_gen_dist <- ggplot(gen_growth_dist_melt, aes(x=GeneticDistance, y=SScoreDiff, col=Condition)) + 
  geom_point() + geom_smooth(method='lm', formula = y~x, level = 0)

#[gen_growth_dist_melt$Condition %in% grep('glucose', unique(gen_growth_dist_melt$Condition), value = TRUE), ]

fits <- lapply(unique(gen_growth_dist_melt$Condition), function(x){
  lm(SScoreDiff ~ GeneticDistance, data = subset(gen_growth_dist_melt, Condition==x))
  })

fit_df <- data.frame(condition = unique(gen_growth_dist_melt$Condition),
                     slope = sapply(fits, function(x){x$coefficients[2]}),
                     int = sapply(fits, function(x){x$coefficients[1]}))

# Lm coefficients are not normally distributed
shapiro.test(fit_df$slope) 
# Shapiro-Wilk normality test
# 
# data:  fit_df$coef
# W = 0.94616, p-value = 0.04324
shapiro.test(fit_df$int) 
# Shapiro-Wilk normality test
# 
# data:  fit_df$int
# W = 0.94483, p-value = 0.03868

# Model variance of scorediff
bins=50
gen_growth_dist_melt$bins <- cut(gen_growth_dist_melt$GeneticDistance, bins, labels = FALSE)

gen_growth_binned <- matrix(data = 1:bins, nrow = bins, ncol = length(unique(gen_growth_dist_melt$Condition)))
colnames(gen_growth_binned) <- unique(gen_growth_dist_melt$Condition)
gen_growth_binned <- melt(gen_growth_binned)[,c(1,2)]
colnames(gen_growth_binned) <- c('Bin', 'Condition')

gen_growth_binned$varSScoreDiff <- NA
for (con in unique(gen_growth_binned$Condition)){
  for (bin in 1:bins){
    gen_growth_binned[gen_growth_binned$Condition == con & gen_growth_binned$Bin == bin, 'varSScoreDiff'] = 
      sd(gen_growth_dist_melt[gen_growth_dist_melt$bins == bin & gen_growth_dist_melt$Condition == con, 'SScoreDiff'])
  }
}
r <- range(gen_growth_dist_melt$GeneticDistance)
gen_growth_binned$Midpoint <- r[1] + gen_growth_binned$Bin * (r[2] - r[1])/(bins * 2)

p_binned_var <- ggplot(gen_growth_binned, aes(x=Midpoint, y=varSScoreDiff, col=Condition)) + geom_point() + 
  geom_smooth(method = 'lm', formula = y~x, level=0) + xlab('Genetic Distance') + ylab("Standard Deviation of S-Score Differences")

fit <- lm(formula = varSScoreDiff ~ Midpoint:Condition, data = gen_growth_binned)

p_binned_sd_box <- ggplot(gen_growth_binned, aes(x=Condition, y=varSScoreDiff)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Standard Deviation of S-Score Differences")

ggsave('figures/condition-growth-diff-sd-box.pdf', device = 'pdf', plot = p_binned_sd_box, width = 12, height = 10)
ggsave('figures/condition-growth-diff-sd.pdf', device = 'pdf', plot = p_binned_var, width = 12, height = 10)

## Test high genetic difference strains only
gen_growth_dist_melt_far <- subset(gen_growth_dist_melt, GeneticDistance > 100000)

p_far_growth_distribution <- ggplot(gen_growth_dist_melt_far, aes(x=Condition, y=SScoreDiff)) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# S-Score normalises to condition
p_growth_box <- ggplot(growth_melt, aes(reorder(Condition, SScore, sd), SScore)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

p_growth_dist_box <- ggplot(gen_growth_dist_melt, aes(reorder(Condition, SScore, sd), SScoreDiff)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Growth differences not all from same population
kruskal.test(SScoreDiff~Condition, data = gen_growth_dist_melt)

# Apply KS-Test over different condition pairs and correct for multiple testing
ks_tests <- sapply(unique(gen_growth_dist_melt$Condition), function(x){
  sapply(unique(gen_growth_dist_melt$Condition), function(y){
    if (!x == y){
      return(ks.test(gen_growth_dist_melt[gen_growth_dist_melt$Condition==x, 'SScoreDiff'],
                     gen_growth_dist_melt[gen_growth_dist_melt$Condition==y, 'SScoreDiff'])$p.value)
    } else {
      return(1)
    }
  })
}) 
dimnames(ks_tests) <- list(unique(gen_growth_dist_melt$Condition),unique(gen_growth_dist_melt$Condition))

ks_tests[lower.tri(ks_tests, diag = TRUE)] <- NA

ks_tests_melt <- subset(melt(ks_tests, value.name = 'p.value', variable.names=c('Cond1', 'Cond2')), !is.na(p.value))
ks_tests_melt$p.adj <- p.adjust(ks_tests_melt$p.value, method = 'bonferroni')

ks_tests_melt <- ks_tests_melt[order(ks_tests_melt$p.adj),]
ks_tests_melt_sig <- subset(ks_tests_melt, p.adj < 0.05)
