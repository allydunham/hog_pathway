setwd("Projects/hog/")

growth <- read.table("data/raw/bede_2017_parsed.tsv",header=TRUE,sep='\t')

strains <- read.table('data/hog-gene-variants.genotypes', header=TRUE, sep='\t', row.names = 1)

hist(colSums(strains))
hist(rowSums(strains), breaks = 2022, xlim = c(0,10))
