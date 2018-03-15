# Script looking at ko probability of genes comparaed to growth rates
setwd('~/Projects/hog/')

ko_prob <- read.table("data/hog-gene-variants.koProb", header = TRUE, row.names = 1, sep='\t')
hog_meta <- read.table("meta/hog-gene-loci")
colnames(hog_meta) <- c("chr","start","end","id","gene","strand")
ko_prob <- ko_prob[,colnames(ko_prob) %in% hog_meta$id]

growth <- read.table("data/raw/bede_2017_parsed.tsv",header=TRUE,sep='\t')

meta <- read.table('meta/strain_information.tsv', sep = '\t', header=TRUE, fill=TRUE, quote = "")
strain_to_sys <- meta$Standardized.name
names(strain_to_sys) <- meta$Isolate.name

growth <- growth[!is.na(strain_to_sys[rownames(growth)]),c("sodium.chloride.0.4mM", "sodium.chloride.0.6mM", "Sorbitol.1mM")]
rownames(growth) <- strain_to_sys[rownames(growth)]

# Gene Correlations
library(gplots)

ko_prob_cor <- ko_prob[,!colSums(ko_prob) == 0]
gene_cor <- cor(as.matrix(ko_prob_cor))
colnames(gene_cor) <- sapply(colnames(gene_cor), function(x){hog_meta[hog_meta$id == x, "gene"]})
rownames(gene_cor) <- sapply(rownames(gene_cor), function(x){hog_meta[hog_meta$id == x, "gene"]})

cols <- colorRampPalette(c("red","white","blue"))(256)
heatmap.2(gene_cor, symm=TRUE, col=cols, breaks=seq(-1,1,2/256), trace = "none")

# Filter to strains with growth data
ko_prob <- ko_prob[rownames(growth),]

ko_prob$growth <- growth$sodium.chloride.0.6mM

for (i in colnames(ko_prob)){
  plot(ko_prob[,i], ko_prob$growth, xlab=i, ylab = "S-Score", pch=20)
  readline(prompt="Press [enter] to continue")
}

fit <- lm(growth ~ ., data = ko_prob)

m <- colMeans(ko_prob)
names(m) <- sapply(names(m), function(x){hog_meta[hog_meta$id == x, "gene"]})
