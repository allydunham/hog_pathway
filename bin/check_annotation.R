# Script to correct problems in vcf annotation
setwd('~/Projects/hog/')

# currently just filter badly annotated genes completely
t <- read.table('data/all-genes-no-missing.impact', header = TRUE, sep='\t')
tt <- t[(duplicated(t[c("mut_id","gene")]) | duplicated(t[c("mut_id","gene")], fromLast = TRUE)), ]

t <- t[!t$gene %in% unique(tt$gene),]

write.table(t, file = 'data/all-genes-no-missing-filtered.impact', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
