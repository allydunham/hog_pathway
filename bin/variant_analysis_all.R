setwd("Projects/hog/")
library(tidyr)

growth <- read.table("data/raw/bede_2017_parsed.tsv",header=TRUE,sep='\t')

strains <- read.table('data/all-sig-genes.genotypes', header=TRUE, sep='\t', row.names = 1)
meta <- read.table('meta/strain_information.tsv', sep = '\t', header=TRUE, fill=TRUE, quote = "")

strain_to_sys <- meta$Standardized.name
names(strain_to_sys) <- meta$Isolate.name

## Filter growth to strains with correct info
growth <- growth[!is.na(strain_to_sys[rownames(growth)]),]
rownames(growth) <- strain_to_sys[rownames(growth)]

## Filter strain genotypes to strains with growth data
strains <- strains[,rownames(growth)]
strains <- strains[!(rowSums(strains) == 0),]

strains.rare <- strains[rowSums(strains) < 50,]

## Summary plots
hist(colSums(strains.rare))
hist(rowSums(strains.rare), breaks = 98)

## Import impact data
impact <- read.table('data/all-sig-genes.impact', sep='\t', header=TRUE)
impact <- impact[impact$mut_id %in% rownames(strains.rare),]

## Normalise impacts relative to ko score
ko_growth <- read.table('data/all-sig-genes.ko', sep='\t', header = TRUE)


## Analyse genotypes
get_con <- function(x){
  return(unique(ko_growth[ko_growth$gene == x, "condition"]))
}

impact$condition <- sapply(impact$gene, get_con)

## Lookup table for condition names
con_hash <- c("2,4-Dichlorophenoxyacetic acid (48H)" = NA,
              "Acetic acid (48H)" = "pH.3", ## Best guess?
              "Amphotericin B (48H)" = NA,
              "Amphotericin B + anaerobic (48H)" = NA,
              "Anaerobic growth (48H)" = "anaerobic",
              "Cadmium chloride (48H)" = NA,
              "Caffeine 15mM (48H)" = "Caffeine.2mM", # Closest match?
              "Caspofungin (72H)" = "caspofungin.100nM",
              "Clozapine (48H)" = NA,
              "Cyclohexamide (48H)" = "Cyclohexamide.100nM",
              "DMSO 1%  (48H)" = "DMSO.1.",
              "Doxorubicin (48H)" = NA,
              "Glucose 20% (48H)" = "High.glucose.20.",
              "Glycerol 2%  (48H)" = "glycerol.2.",
              "Glycerol 2%  (72H)" = "glycerol.2.",
              "Maltose 2%  (48H)" = NA,
              "Maltose 2%  (72H)" = NA,
              "NaCl 0.4M (48H)" = "sodium.chloride.0.4mM",
              "NaCl 0.4M + 39ºC (48H)" = NA,
              "NaCl 0.4M + 39ºC (72H)" = NA,
              "NaCl 0.6M (48H)" = "sodium.chloride.0.6mM",
              "NaCl 0.6M (72H)" = "sodium.chloride.0.6mM",
              "NaCl 0.6M + 39ºC (48H)" = NA,
              "NaCl 0.6M + 39ºC (72H)" = NA,
              "NiSO4 (48H)" = NA,
              "Nitrogen starvation (48H)" = "nitrogen.starvation",
              "Nystatin (48H)" = "nystatin.30uM",
              "Paraquat (48H)" = "paraquat.1.2mM",
              "Paraquat (72H)" = "paraquat.1.2mM",
              "SC + hepes (48H)" = "Synthetic.complete.media",
              "Sorbitol 1M (48H)" = "Sorbitol.1mM",
              "aa starvation (48H)" = "amino.acid.depravation..HIS.LEU.URA.LYS.")

impact <- impact[sapply(impact$condition, length) == 1,] # Only interest in effect that impact a single condition for simplicity
impact$condition <- unlist(impact$condition)

impact <- impact[!is.na(con_hash[impact$condition]),]

get_growth <- function(x){
  s <- colnames(strains)[which(strains[x["mut_id"],] > 0)]
  t <- growth[s, con_hash[x["condition"]]]
  names(t) <- s
  return(t)
}

impact$growth <- apply(impact, 1, get_growth)
impact <- unnest(impact, growth)

impact$norm_growth <- impact$growth/unlist(apply(impact, 1, function(x){
  mean(ko_growth[ko_growth$condition == unlist(x["condition"]) & ko_growth$gene == unlist(x["gene"]),"score"])
}))

hist(impact$sift_score)
plot(impact$sift_score, impact$norm_growth) # No apparant relation with sift score
boxplot(norm_growth ~ type, data=impact) # Possible relation with frameshift and nonsense
plot(impact$pos_aa, impact$norm_growth) # No apparant relation with sift score
boxplot(norm_growth ~ as.factor(ref_codon), data=impact)
boxplot(norm_growth ~ as.factor(alt_aa), data=impact[nchar(impact$alt_aa) == 1,])
plot(impact$foldx_ddG, impact$norm_growth) # Nothing strong found on any factor

fit <- lm(norm_growth ~ sift_score, data=impact) ## Sift score is significant but irrelevant

### Test ML solutions
library(caret)
impact.ml <- impact[!is.na(impact$sift_score) & !is.na(impact$norm_growth),]

train_index <- createDataPartition(impact.ml$growth, times = 1, p = 0.8, list = FALSE)
training <- impact.ml[train_index,]
test <- impact.ml[-train_index,]

model <- train(norm_growth ~ type + sift_score + ref_codon + alt_codon + ref_aa + alt_aa,
               data=training, method = "bayesglm")

prediction <- predict(model, newdata = test)

plot(test$norm_growth, prediction)

