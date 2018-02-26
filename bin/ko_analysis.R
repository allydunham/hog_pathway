setwd("Projects/pathways/variation/")

## Import gene ko growth scores for hog genes and compare impacts
hog_ko_growth <- read.table("hog_ko_scores.tsv",sep = "\t")
colnames(hog_ko_growth) <- c("strain","condition","geneID","geneName","position","score","qval")

## Limit to standard lab strain
hog_ko_growth <- hog_ko_growth[hog_ko_growth$strain == "S288C",]
hog_ko_growth$col <- "black"
hog_ko_growth[hog_ko_growth$qval < 0.05,"col"] <- "red"

osmotic_shock <- grep("NaCl|Sorbitol",unique(hog_ko_growth$condition),value = TRUE)

for (con in osmotic_shock){
  t <- hog_ko_growth[hog_ko_growth$condition == con,]
  plot(NA, main = con, ylim = c(0,1),xlim = c(-8,8))
  text(t$score,t$qval,col=t$col,labels = t$geneName,cex = 0.75)
  readline(prompt="Press [enter] to continue")
}

## Same analysis with all genes
ko_growth <- read.table("osmotic_ko_scores.tsv",sep = "\t",quote="")
colnames(ko_growth) <- c("strain","condition","geneID","geneName","position","score","qval")
ko_growth <- ko_growth[!ko_growth$geneID == "WT",]

ko_growth <- ko_growth[ko_growth$strain == "S288C",]
ko_growth$col <- "black"
ko_growth[ko_growth$qval < 0.05,"col"] <- "red"

osmotic_shock <- unique(ko_growth$condition)
osmotic_shock <- osmotic_shock[c(1,4,5,8)]
for (con in osmotic_shock){
  t <- ko_growth[ko_growth$condition == con,]
  plot(NA, main = con, ylim = c(0,1),xlim = c(-8,8))
  text(x = t$score, y = t$qval, labels = t$geneName, col=t$col,cex = 0.75)
  readline(prompt="Press [enter] to continue")
}

sig <- ko_growth[ko_growth$qval < 0.05,]
sigl <- unlist(sapply(osmotic_shock, function(x){unique(sig[sig$condition == x,"geneID"])}))
table(sigl)
# cat(unique(sigl),sep = "\n") # <- for GO enrichment
