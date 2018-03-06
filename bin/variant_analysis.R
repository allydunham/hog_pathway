setwd("Projects/hog/")
library(tidyr)

growth <- read.table("data/raw/bede_2017_parsed.tsv",header=TRUE,sep='\t')

strains <- read.table('data/hog-gene-variants.genotypes', header=TRUE, sep='\t', row.names = 1)
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
impact <- read.table('data/hog-gene-variants.impact', sep='\t', header=TRUE)
impact <- impact[impact$mut_id %in% rownames(strains.rare),]



## Analyse genotypes
get_growth <- function(id, con){
  s <- colnames(strains)[which(strains[id,] > 0)]
  return(growth[s, con])
}

impact$growth <- sapply(impact$mut_id, get_growth, con="sodium.chloride.0.6mM")
impact <- unnest(impact, growth)

## Normalise impacts relative to ko score
ko_growth <- read.table('data/hog_ko_scores.tsv', sep='\t')
colnames(ko_growth) <- c('strain', 'condition', 'id', 'gene', 'pos', 's-score', 'p-val')
ko_growth <- ko_growth[ko_growth$strain == "S288C" & ko_growth$condition == "NaCl 0.6M (48H)",]

impact$norm_growth <- impact$growth/sapply(impact$gene, function(x){mean(ko_growth[ko_growth$id == x, "s-score"])})

hist(impact$sift_score)
plot(impact$sift_score, impact$norm_growth) # No apparant relation with sift score
boxplot(norm_growth ~ type, data=impact) # Possible relation with frameshift and nonsense
boxplot(norm_growth ~ as.factor(ptm_modification), data=impact)
plot(impact$foldx_ddG, impact$norm_growth, ylim=c(-20,20)) # Nothing strong found on any factor


### Test ML solutions
library(caret)
impact.ml <- impact[!is.na(impact$sift_score) & !is.na(impact$norm_growth),]

train_index <- createDataPartition(impact.ml$growth, times = 1, p = 0.8, list = FALSE)
training <- impact.ml[train_index,]
test <- impact.ml[-train_index,]

model <- train(norm_growth ~ type + sift_score + ref_codon + alt_codon + ref_aa + alt_aa,
               data=training, method = "rf")

prediction <- predict(model, newdata = test)



