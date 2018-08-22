# Script examining bede growth data
setwd('Projects/pathways/variation/')

growth <- read.table("raw_data/bede_2017_parsed.tsv",sep='\t',header = TRUE)
growth_reps <- read.table("raw_data/bede_2017_parsed_rep1_rep2.tsv",sep='\t',header = TRUE)

growth_rep1 <- growth_reps[,1:43]
growth_rep2 <- growth_reps[,44:86]

means <- sapply(1:43, function(con){sapply(1:166,function(strain){
  mean(c(growth_rep1[strain,con],growth_rep2[strain,con]))
})})
rownames(means) <- rownames(growth)
colnames(means) <- colnames(growth)

t <- means == growth
166*43 - sum(t)

plot(unname(unlist(growth_rep1)),unname(unlist(growth_rep2)),pch=20,
     main = "Growth Rate Reps",xlab = "Rep1",ylab = "Rep2")
abline(0,1)


diff <- unlist(growth_rep1 - growth_rep2)
hist(diff,xlim = c(-4,4),breaks = 30,main = "Histogram of difference between reps",xlab = "Difference")
