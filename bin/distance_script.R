## Script to calculate the genetic distance between strains and output comparisons to growth

args = commandArgs(trailingOnly = TRUE)

genotype_file = args[1] # matrix of genotypes (0/1/2) of each variant allele
growth_file = args[2] # Table of growth S-scores per condition
meta_file = args[3] # Table of strain meta info

# Import growth file
growth <- read.table(growth_file,header=TRUE,sep='\t')
growth <- growth[,c("sodium.chloride.0.4mM", "sodium.chloride.0.6mM")]

# Import meta file giving translation between common (used for growth data) and standardised name (used in vcf)
meta <- read.table(meta_file, sep = '\t', header=TRUE, fill=TRUE, quote = "")
strain_to_sys <- meta$Standardized.name
names(strain_to_sys) <- meta$Isolate.name

growth <- growth[!is.na(strain_to_sys[rownames(growth)]),]
rownames(growth) <- strain_to_sys[rownames(growth)]

# Import genotypes
genotypes <- t(read.table(genotype_file, header=TRUE, row.names = 1))
genotypes_growth <- genotypes[rownames(growth),] # Filter to those strains with growth data

# Calculate genetic distance
genetic_distance <- as.matrix(dist(genotypes, method = "man"))

genetic_distance_growth <- as.matrix(dist(genotypes_growth, method = "man"))
growth_distance_NaCl4 <- as.matrix(dist(growth[,1]))
growth_distance_NaCl6 <- as.matrix(dist(growth[,2]))

# Save copy of calculated distance data
save(genetic_distance, genetic_distance_growth, growth_distance_NaCl4, growth_distance_NaCl6, file = 'genetic_distance.Rdata')

# Save diagnostic plots
pdf('genetic_dist_nacl0.4mM.pdf', 12, 8)
plot(genetic_distance_growth[upper.tri(genetic_distance)], growth_distance_NaCl4[upper.tri(growth_distance_NaCl4)],
     pch=20, xlab="Genetic Distance (Manhatten)", ylab="Difference in S-Score (NaCl 0.4mM)",
     ylim=c(0,7), main='Effect of genetic distance on growth rate in 0.4mM NaCl')
dev.off()

pdf('genetic_dist_nacl0.6mM.pdf', 12, 8)
plot(genetic_distance_growth[upper.tri(genetic_distance)], growth_distance_NaCl6[upper.tri(growth_distance_NaCl6)],
     pch=20, xlab="Genetic Distance (Manhatten)", ylab="Difference in S-Score (NaCl 0.6mM)",
     ylim=c(0,7), main='Effect of genetic distance on growth rate in 0.6mM NaCl')
dev.off()

