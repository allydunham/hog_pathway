setwd("~/Projects/hog/")
library(Biostrings)

genotypes <- t(read.table("data/hog-gene-variants.all-genotypes", header=TRUE, row.names = 1))

growth <- read.table("data/raw/bede_2017_parsed.tsv",header=TRUE,sep='\t')
growth <- growth[,c("sodium.chloride.0.4mM", "sodium.chloride.0.6mM", "Sorbitol.1mM")]

meta <- read.table('meta/strain_information.tsv', sep = '\t', header=TRUE, fill=TRUE, quote = "", comment.char = '')
strain_to_sys <- meta$Standardized.name
names(strain_to_sys) <- meta$Isolate.name
sys_to_strain <- names(strain_to_sys)
names(sys_to_strain) <- strain_to_sys

growth <- growth[!is.na(strain_to_sys[rownames(growth)]),]
rownames(growth) <- strain_to_sys[rownames(growth)]

genotypes <- genotypes[rownames(growth),]

genetic_distance <- as.matrix(dist(genotypes, method = "man"))
growth_distance_NaCl4 <- as.matrix(dist(growth[,1]))
growth_distance_NaCl6 <- as.matrix(dist(growth[,2]))
growth_distance_Sorb <- as.matrix(dist(growth[,3]))

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
growth_distance <- lapply(colnames(growth), function(x){as.matrix(dist(growth[,x]))})
names(growth_distance) <- colnames(growth)

for (i in 1:length(growth_distance)){
  colnames(growth_distance[[i]]) <- rownames(growth)
  rownames(growth_distance[[i]]) <- rownames(growth)
}
