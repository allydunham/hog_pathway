# Script looking at ko probability of genes comparaed to growth rates
setwd('~/Projects/hog/')

ko_prob <- read.table("data/hog-gene-variants.probs.blosum", header = TRUE, row.names = 1, sep='\t')
hog_meta <- read.table("meta/hog-gene-loci")
colnames(hog_meta) <- c("chr","start","end","id","gene","strand")
ko_prob <- ko_prob[,colnames(ko_prob) %in% hog_meta$id]
colnames(ko_prob) <- sapply(colnames(ko_prob), function(x){hog_meta[hog_meta$id == x, "gene"]})

growth <- read.table("data/raw/bede_2017_parsed.tsv",header=TRUE,sep='\t')
rownames(growth) <- toupper(rownames(growth))

meta <- read.table('meta/strain_information.tsv', sep = '\t', header=TRUE, fill=TRUE, quote = "", comment.char = "")
meta$Isolate.name <- toupper(meta$Isolate.name)
strain_to_sys <- meta$Standardized.name
names(strain_to_sys) <- meta$Isolate.name

growth <- growth[!is.na(strain_to_sys[rownames(growth)]),c("sodium.chloride.0.4mM", "sodium.chloride.0.6mM", "Sorbitol.1mM")]
rownames(growth) <- strain_to_sys[rownames(growth)]

pdf('figures/p_affected_genes_box_corrected.pdf', width = 12, height = 10)
boxplot(ko_prob, ylab="P(Affected)", las=2, main='Per gene impact probability distributions', pch=20)
dev.off()

# Gene Correlations
library(gplots)

ko_prob_cor <- ko_prob[,!colSums(ko_prob) == 0]
gene_cor <- cor(as.matrix(ko_prob_cor))
#colnames(gene_cor) <- sapply(colnames(gene_cor), function(x){hog_meta[hog_meta$id == x, "gene"]})
#rownames(gene_cor) <- sapply(rownames(gene_cor), function(x){hog_meta[hog_meta$id == x, "gene"]})

cols <- colorRampPalette(c("red","white","blue"))(256)

pdf('figures/gene_ko_probs_heatmap_corrected.pdf', width = 12, height = 12)
heatmap.2(gene_cor, symm=TRUE, col=cols, breaks=seq(-1,1,2/256), trace = "none")
dev.off()

# Filter to strains with growth data
ko_prob <- ko_prob[rownames(growth),]

ko_prob$growth <- growth$sodium.chloride.0.6mM

for (i in colnames(ko_prob)){
  plot(ko_prob[,i], ko_prob$growth, xlab="P(Affected)", main=i, ylab = "S-Score", pch=20)
  readline(prompt="Press [enter] to continue")
}

fit <- lm(growth ~ ., data = ko_prob)

m <- colMeans(ko_prob)
names(m) <- sapply(names(m), function(x){hog_meta[hog_meta$id == x, "gene"]})

### Analyse whole pathway ko prob
path_active <- read.table('data/hog-gene-variants.path.blosum', header = TRUE, sep='\t')
path_active <- subset(path_active, strain %in% rownames(growth))
rownames(path_active) <- path_active$strain
path_active <- merge(path_active, growth, by='row.names')
path_active$Row.names <- NULL

pdf('figures/bin-path-vs-growth.pdf', width = 12, height = 10)
boxplot(sodium.chloride.0.6mM ~ hog_active, data = path_active, ylab = 'S-Score', xlab = 'HOG Path Active',
        main = 'Binary pathway activty impact on growth')
dev.off()

plot(path_active$hog_active, path_active$sodium.chloride.0.4mM)

t.test(sodium.chloride.0.6mM ~ hog_active, data = path_active, alternative='less')

## Hog1 activity probability against growth
pdf('figures/path-probability-vs-growth.pdf', width = 12, height = 10)
plot(path_active$hog_probability, path_active$sodium.chloride.0.4mM, pch=20,
     xlab = 'Hog1 Activity Probability', ylab = 'S-Score', main = 'Effect of pathway activity on growth')
fit4 <- lm(sodium.chloride.0.4mM ~ hog_probability, data = path_active)
abline(a = fit4$coefficients[1], b = fit4$coefficients[2])

points(path_active$hog_probability, path_active$sodium.chloride.0.6mM, pch=20, col='red')
fit6 <- lm(sodium.chloride.0.6mM ~ hog_probability, data = path_active)
abline(a = fit6$coefficients[1], b = fit6$coefficients[2], col='red')

legend('bottomright', legend = c('0.4mM NaCl', '0.6mM NaCl'), fill = c('black', 'red'), bty = 'n')
dev.off()

## growth against growth
pdf('figures/growth_correlation_nacl.pdf', width = 12, height = 10)
plot(path_active$sodium.chloride.0.4mM, path_active$sodium.chloride.0.6mM, pch=20, xlim = c(-4,4), ylim=c(-4,4),
     xlab = 'S-Score (NaCl 0.4mM)', ylab = 'S-Score (NaCl 0.6mM)', main = 'Growth Correlation')
fitG <- lm(sodium.chloride.0.4mM ~ sodium.chloride.0.6mM, data = path_active)
abline(a = fitG$coefficients[1], b = fitG$coefficients[2])
dev.off()


