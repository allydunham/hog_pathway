setwd("~/Projects/hog/")
library(Biostrings)

genotypes <- t(read.table("data/hog-gene-variants.genotypes", header=TRUE, row.names = 1))

growth <- read.table("data/raw/bede_2017_parsed.tsv",header=TRUE,sep='\t')
growth <- growth[,c("sodium.chloride.0.4mM", "sodium.chloride.0.6mM", "Sorbitol.1mM")]

meta <- read.table('meta/strain_information.tsv', sep = '\t', header=TRUE, fill=TRUE, quote = "")
strain_to_sys <- meta$Standardized.name
names(strain_to_sys) <- meta$Isolate.name

growth <- growth[!is.na(strain_to_sys[rownames(growth)]),]
rownames(growth) <- strain_to_sys[rownames(growth)]

genotypes <- genotypes[rownames(growth),]

genetic_distance <- as.matrix(dist(genotypes, method = "man"))
growth_distance_NaCl4 <- as.matrix(dist(growth[,1]))
growth_distance_NaCl6 <- as.matrix(dist(growth[,2]))
growth_distance_Sorb <- as.matrix(dist(growth[,3]))

plot(genetic_distance[upper.tri(genetic_distance)], growth_distance_NaCl4[upper.tri(growth_distance_NaCl4)],
     pch=20, xlab="Genetic Distance (Manhatten)", ylab="Difference in S-Score (NaCl 0.4mM)")

plot(genetic_distance[upper.tri(genetic_distance)], growth_distance_NaCl6[upper.tri(growth_distance_NaCl6)],
     pch=20, xlab="Genetic Distance (Manhatten)", ylab="Difference in S-Score (NaCl 0.6mM)")

plot(genetic_distance[upper.tri(genetic_distance)], growth_distance_Sorb[upper.tri(growth_distance_Sorb)],
     pch=20, xlab="Genetic Distance (Manhatten)", ylab="Difference in S-Score (Sorbitol 1mM)")








