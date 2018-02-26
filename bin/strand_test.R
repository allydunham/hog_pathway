t <- read.table("data/hog-gene-variants.genomic",header = TRUE)
t <- t$mut_id
tt <- read.table("data/hog-gene-variants.function", header=TRUE, sep='\t')

ttt <- t[!t %in% coding$mut_id]
