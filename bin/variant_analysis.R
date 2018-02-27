setwd("Projects/hog/")

growth <- read.table("data/raw/bede_2017_parsed.tsv",header=TRUE,sep='\t')

strains <- read.table('data/combined-gene-variants.genotypes', header=TRUE, sep='\t', row.names = 1)
meta <- read.table('meta/strain_information.tsv', sep = '\t', header=TRUE, fill=TRUE, quote = "")

strain_to_sys <- meta$Standardized.name
names(strain_to_sys) <- meta$Isolate.name

## Filter growth to strains with correct info
growth <- growth[!is.na(strain_to_sys[rownames(growth)]),c("sodium.chloride.0.4mM", "sodium.chloride.0.6mM", "Sorbitol.1mM")]
rownames(growth) <- strain_to_sys[rownames(growth)]

## Filter strain genotypes to strains with growth data
strains <- strains[,rownames(growth)]
strains <- strains[!(rowSums(strains) == 0),]

strains.rare <- strains[rowSums(strains) < 50,]

## Summary plots
hist(colSums(strains.rare))
hist(rowSums(strains.rare), breaks = 98)

## Import impact data
impact <- read.table('data/combined-gene-variants.impact', sep='\t', header=TRUE)
impact <- impact[impact$mut_id %in% rownames(strains.rare),]

## Analyse genotypes
get_growth <- function(id, con){
  s <- colnames(strains)[which(strains[id,] > 0)]
  return(mean(growth[s, con]))
}

impact$growth <- sapply(impact$mut_id, get_growth, con="sodium.chloride.0.6mM")

hist(impact$sift_score)
plot(impact$sift_score, impact$growth) # No apparant relation with sift score
boxplot(growth ~ type, data=impact) # Possible relation with frameshift and nonsense
boxplot(growth ~ as.factor(ptm_modification), data=impact)
plot(impact$foldx_ddG, impact$growth) # Nothing strong found on any factor



