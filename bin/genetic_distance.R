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

genetic_distance <- as.matrix(dist(genotypes, method = "manhattan"))
growth_distance <- as.matrix(dist(growth[,1]))

plot(genetic_distance[upper.tri(genetic_distance)], growth_distance[upper.tri(growth_distance)])


## Genotype against growth
impact <- read.table('data/hog-gene-variants.impact', sep='\t', header=TRUE)
impact <- impact[impact$mut_id %in% colnames(genotypes),]

# Only consider non-synonymous coding changes (need to reduce number of variables to consider)
impact <- impact[!impact$type == "synonymous",]
impact <- impact[!is.na(impact$id),]

data("BLOSUM62")
impact$blosum <- apply(impact, 1, function(x){
  if (x["alt_aa"] %in% rownames(BLOSUM62)){
    return(BLOSUM62[unlist(x["ref_aa"]),unlist(x["alt_aa"])])
  } else {
    return(NA)
  }
  })

## Only consider subs that generally occur less than expected
impact <- impact[is.na(impact$blosum) | impact$blosum < 0,]

genotypes <- genotypes[,unique(impact$mut_id)]
genotypes <- genotypes[,colSums(genotypes) < dim(genotypes)[1]*2*0.1]

gen_df <- as.data.frame(genotypes)
gen_df$growth <- as.numeric(growth$sodium.chloride.0.6mM)

library(caret)
train_index <- createDataPartition(gen_df$growth, times = 1, p = 0.8, list = FALSE)
training <- gen_df[train_index,]
test <- gen_df[-train_index,]
  
model <- train(growth ~ ., data = training, method = "rf")

model <- train(growth ~ ., data = training, method = "bayesglm")


prediction <- predict(model, newdata = test)
plot(test$growth, prediction)













