## Script exploring the distribution of neutral probabilities based on data from Jelier et al. 2011
setwd('~/Projects/hog/')
library(Biostrings)
data("BLOSUM62")

impact <- read.table("data/jelier.impact", sep = '\t', header=TRUE)
impact$blosum <- apply(impact, 1, function(x){
  if (x["alt_aa"] %in% rownames(BLOSUM62) & x["ref_aa"] %in% rownames(BLOSUM62)){
    return(BLOSUM62[unlist(x["ref_aa"]),unlist(x["alt_aa"])])
  } else if (x["alt_aa"] == 'O' | x["ref_aa"] == 'O'){
    return(-4) # Set any with an unusual AA as impactful, under the assumption any use of these is meaningful
  } else{
    return(NA)
  }
})

impact$type <- NA
impact[!impact$alt_aa == impact$ref_aa, 'type'] <- 'nonsynonymous'
impact[impact$alt_aa == impact$ref_aa, 'type'] <- 'synonymous' # Non synonymous and only one nonsense mutation in training set, maybe doesn't matter as could be dealt with separately?
impact[impact$alt_aa == '*', 'type'] <- 'nonsense'
impact$effect <- as.factor(impact$effect)

boxplot(sift_score ~ effect, data = impact) # Impactful mutations clearly have much lower sift score
plot(impact$foldx_ddG, impact$effect)

boxplot(foldx_ddG ~ effect, data = impact)

library(caret)
library(e1071)
t <- impact[!is.na(impact$sift_score),]
train_index <- createDataPartition(t$effect, times = 1, p = 0.8, list = FALSE)
training <- t[train_index,]
test <- t[-train_index,]

model <- train(effect ~ pos_aa + blosum + sift_score + elm_lost, data = training, method = "rf")
## using blosum gives small improvement over ref/alt aa

### Full version with processing - no meaningfull improvement
preProc <- preProcess(t, method = c("center","scale"))
t_pro <- predict(preProc, t)

train_index <- createDataPartition(t_pro$effect, times = 1, p = 0.8, list = FALSE)
training <- t_pro[train_index,]
test <- t_pro[-train_index,]

tc <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10)

model <- train(effect ~ pos_aa + blosum + sift_score + elm_lost, data = training, method = "rf")

prediction <- predict(model, newdata = test)
table(test$effect, prediction)
summary(model)
varImp(object=model)

probs <- predict(model, newdata = test, type = "prob")
boxplot(probs[,2] ~ test$effect)

## Ranger implementation of RF
library(ranger)

model = ranger(effect ~ pos_aa + blosum + sift_score + elm_lost, data = training, probability = TRUE)
probs <- predict(model, data = test)
boxplot(probs$predictions[,2] ~ test$effect)

### Fit model of sift to prob instead
t <- impact[!is.na(impact$sift_score),]
t$bin <- cut(t$sift_score, breaks = 20, labels = FALSE)
p <- sapply(1:20, function(x){
  sum(t[t$bin == x, "effect"] == "0")/sum(t$bin == x)
})

t$p <- p[t$bin]

min_sift <- min(t[t$sift_score > 0,"sift_score"])

f <- function(x,a,b){1/(1 + exp(a*log(x + min_sift) + b))}
fit_sift <- nls(p ~ f(sift_score,a,b), data = t, start = list(a=-1,b=1))

plot(t$sift_score, t$p, pch=20, main = "Probability of Neutrality (SIFT)", xlab = "Sift Score", ylab = "p", ylim=c(0,1))
curve(f(x,coef(fit_sift)[1],coef(fit_sift)[2]),add = TRUE)

# Same workflow on foldx
t <- impact[!is.na(impact$foldx_ddG),]
t$bin <- cut(t$foldx_ddG, breaks = 20, labels = FALSE)
p <- sapply(1:20, function(x){
  sum(t[t$bin == x, "effect"] == "0")/sum(t$bin == x)
})

t$p <- p[t$bin]
f <- function(x,a,b){1/(1 + exp(a*x + b))}
fit_foldx <- nls(p ~ f(foldx_ddG,a,b), data = t, start = list(a=0, b=0))

plot(t$foldx_ddG, t$p, ylim = c(0,1), pch=20, main = "Probability of Neutrality (FoldX)", xlab = "ddG", ylab = "p")
curve(f(x,coef(fit_foldx)[1],coef(fit_foldx)[2]),add = TRUE)

# Same workflow on blosum
t <- impact[!is.na(impact$blosum) & impact$blosum < 5,]
# so few examples of higher blosum scores and all 3 are synonymous changes that must have effects through other mechanisms
t$p <- sapply(t$blosum, function(x){
  sum(t[t$blosum == x, "effect"] == "0")/sum(t$blosum == x)
})

fit_blosum <- lm(p ~ blosum, data = t)

plot(t$blosum, t$p, ylim = c(0,1), pch=20, main = "Probability of Neutrality (Blosum62)", xlab = "Blosum62 Score", ylab = "p")
abline(fit_blosum$coefficients[1], fit_blosum$coefficients[2])
