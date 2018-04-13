# Script documenting filtering applied to ko library to determine significant genes for osmotic shock
setwd('Projects/hog/')
ko <- read.table("data/raw/ko_scores.txt", header=TRUE, sep = '\t', quote = '')
ko <- ko[ko$strain == "S288C",]
ko <- ko[ko$qvalue < 0.05,]
ko <- ko[!ko$gene == "WT",]
write.table(ko, file = "data/all-sig-genes.ko", sep='\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
