# Script to correct problems in vcf annotation
setwd('~/Projects/hog/')

# currently just filter badly annotated genes completely
imp <- read.table('data/all-genes-no-missing.impact', header = TRUE, sep='\t')
dupes <- imp[(duplicated(imp[c("mut_id","gene")]) | duplicated(imp[c("mut_id","gene")], fromLast = TRUE)), ]

imp <- imp[!imp$gene %in% unique(dupes$gene),]

write.table(imp, file = 'data/all-genes-no-missing-filtered.impact', sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
