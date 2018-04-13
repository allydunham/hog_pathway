# Copyright (c) 2017  Genome  Research  Ltd.
# Author: Alistair Dunham
# This  program  is free  software: you  can  redistribute  it and/or  modify  it  under
# the  terms  of the  GNU  General  Public  License  as  published  by the  Free  Software
# Foundation; either  version 3 of the  License , or (at your  option) any  later
# version.
# This  program  is  distributed  in the  hope  that it will be useful , but  WITHOUT
# ANY  WARRANTY; without  even  the  implied  warranty  of  MERCHANTABILITY  or  FITNESS
# FOR A PARTICULAR  PURPOSE. See  the  GNU  General  Public  License  for  more
# details.
# You  should  have  received a copy of the  GNU  General  Public  License  along  with
# this  program. If not , see <http :// www.gnu.org/licenses/>.

## Script containing functions to support QC analysis of 10X variant calls

############################################ Import Librarys and Definitions ########################################
library(e1071)
library(randomForest)
library(fastAdaboost)
library(klaR)
library(caret)
library(data.table)
library(gbm)
library(plyr)

assesmentMetrics <- c("Accuracy","FDR","FNR","FPR","Kappa","Specificity","Sensitivity","Precision","F1")


############################################ Data Import ########################################
loadDels <- function(dat,cla=NA,drop=TRUE){
  dels <- read.table(dat,header=TRUE,
                     colClasses = c("character","factor","integer","integer","integer","factor","numeric","factor","factor",
                                    "integer","integer","numeric","integer","integer","integer","integer","integer","integer","integer"
                                    ,"character","character"))
  
  ## Process read depth
  dels$ReadDepth <- strsplit(as.character(dels$ReadDepth),",")
  dels$ReadDepth <- lapply(dels$ReadDepth,as.numeric)
  
  #poly(seq(0,1,length.out = x$Length + 200),2)
  fits <- apply(dels,1,function(x){w <- seq(0,1,length.out = x$Length + 200) - 0.5 ;lm(x$ReadDepth ~ I(w^2))})
  dels$int <- sapply(fits,function(x){x$coefficients[1]})
  dels$coef1 <- sapply(fits,function(x){x$coefficients[2]})
  #  dels$coef2 <- sapply(fits,function(x){x$coefficients[3]})
  
  dels <- dels[,c("ID","Chromosome","Start","Stop","Length","Ref","Quality","Genotype","CentromereDist","TelomereDist","GC",
                  "LINEs","SINEs","LTRs","LowComplexityRepeats","SimpleRepeats","OtherRepeats","SDs","Source","int","coef1","ReadDepth")]
  
  dels$Source <- as.factor(dels$Source)
  dels$CentBins <- as.factor(dels$CentromereDist < 0.5 * 10^6)
  dels$TeloBins <- as.factor(dels$TelomereDist < 0.5 * 10^6)
  dels$SDbins <- as.factor(dels$SDs > 0)
  
  ## Add annotations
  if (is.na(cla)){
    return(dels)
  } else {
    annot <- read.table(cla,header=TRUE)
    mer <- merge(dels,annot[,c("chr","vcf_start","Manual_check")],by.x = c("Chromosome","Start"),by.y=c("chr","vcf_start"),all.x = TRUE)
    
    if (drop){
      mer <- mer[!is.na(mer$Manual_check),]
    }
    return(mer)
  }
}

############################################ Ensemble Classifier ########################################
## Models is a named list containing trained classifiers to combined into the ensemble
ensembleclassifier <- function(models,train){
  m <- lapply(models,function(x){x(train)})
  class(m) <- c("ensembleClassifier",class(m))
  return(m)
}

## Predictions made by simple majority vote
predict.ensembleClassifier <- function(ens,data){
  p <- apply(sapply(ens,function(x){predict(x,data)}),2,as.numeric)
  return(rowSums(p) > length(ens)/2)
}

############################################ Classifier Tests ########################################
classifier.cohensKappa <- function(m){
  pe <- sum(sapply(1:dim(m)[1],function(x){sum(m[x,])*sum(m[,x])}))/sum(m)^2
  po <- sum(diag(m))/sum(m)
  return((po - pe)/(1 - pe))
}

classifier.accuracy <- function(m){
  return(sum(diag(m))/sum(m))
}

classifier.fdr <- function(m){
  return(m[2,1]/(m[2,1] + m[2,2]))
}

classifier.fnr <- function(m){
  return(m[1,2]/max(m[1,2] + m[2,2],1))
}

classifier.fpr <- function(m){
  return(m[2,1]/max(m[1,1] + m[2,1],1))
}

classifier.precision <- function(m){
  return(m[2,2]/sum(m[2,]))
}

classifier.specificity <- function(m){
  return(m[1,1]/sum(m[,1]))
}

classifier.sensitivity <- function(m){
  return(m[2,2]/sum(m[,2]))
}

classifier.F1 <- function(m){
  tpr <- classifier.sensitivity(m)
  pre <- classifier.precision(m)
  return(2*tpr*pre/(tpr + pre))
}

## Get all metrics about a classifier
classifier.test <- function(model,testData,cat = 'effect'){
  p <- predict(model,testData)
  m <- table(pred=p,true=testData[,cat])
  
  t <- c(classifier.accuracy(m),
         classifier.fdr(m),
         classifier.fnr(m),
         classifier.fpr(m),
         classifier.cohensKappa(m),
         classifier.specificity(m),
         classifier.sensitivity(m),
         classifier.precision(m),
         classifier.F1(m))
  names(t) <- assesmentMetrics
  return(t)
}

## Function to perfrom k-fold cross validation. trainer function must perform training on a single required argument, the training set
classifier.crossValidate <- function(trainer,testData,divs=10){
  ids <- sample.int(dim(testData)[1],dim(testData)[1],FALSE)
  l <- floor(length(ids)/divs)
  r <- length(ids) %% l
  
  g <- as.list(rep(NA,divs))
  c <- 1
  for (i in 1:divs){
    if (i <= r){
      g[[i]] <- ids[c:(c+l)]
      c <- c + l + 1
    } else {
      g[[i]] <- ids[c:(c+l-1)]
      c <- c + l
    }
  }
  
  return(sapply(g,function(x){classifier.test(trainer(testData[-x,]),testData[x,])}))
}

## Identify parameters using forward selection. trainer function takes training data as its argument, using all columns to train
classifier.forwardSelection <- function(trainer,testData,vars,divs=10,metric='Kappa',sel=which.max){
  v <- c()
  use <- rep(TRUE,length(vars))
  names(use) <- vars
  
  resVars <- as.list(rep(NA,length(vars)))
  resMetric <- matrix(NA,nrow = length(vars),ncol = 9)
  colnames(resMetric) <- c("Accuracy","FDR","FNR","FPR","Kappa","Specificity","Sensitivity","Precision","F1")
  
  ## Consider sequentially longer variable sets, adding best outcome each time
  for (i in 1:length(vars)){
    t <- sapply(vars[use],function(x){rowMeans(classifier.crossValidate(trainer,testData[,c("cat",v,x),with=FALSE],divs))})
    n <- vars[use][sel(t[metric,])]
    v <- c(v,n)
    use[n] <- FALSE
    resVars[[i]] <- v
    resMetric[i,] <- t[,n]
  }
  return(list(vars=resVars,metric=resMetric))
}

## Identify parameters using reverse selection. trainer function takes training data as its argument, using all columns to train
classifier.reverseSelection <- function(trainer,testData,vars,divs=10,metric='Kappa',sel=which.max){
  v <- vars
  
  resVars <- as.list(rep(NA,length(vars)))
  resMetric <- matrix(NA,nrow = length(vars),ncol = 9)
  colnames(resMetric) <- c("Accuracy","FDR","FNR","FPR","Kappa","Specificity","Sensitivity","Precision","F1")
  
  resVars[[length(vars)]] <- v
  resMetric[length(vars),] <- rowMeans(classifier.crossValidate(trainer,testData[,c("cat",v),with=FALSE],divs))
  
  ## Consider sequentially longer variable sets, adding best outcome each time
  for (i in (length(vars)-1):1){
    t <- sapply(1:length(v),function(x){rowMeans(classifier.crossValidate(trainer,testData[,c("cat",v[-x]),with=FALSE],divs))})
    n <- sel(t[metric,])
    v <- v[-n]
    resVars[[i]] <- v
    resMetric[i,] <- t[,n]
  }
  return(list(vars=resVars,metric=resMetric))
}

############################################ Classifier Plotting Functions ########################################
## Plot parameter selection results
plotFactors <- function(x,main="",ylim=c(0,1),xlab="Number of Factors",ylab="Factor Value",
                        cols=c("black","firebrick2","cornflowerblue","purple","green","orange","yellow","hotpink","lightblue")){
  par(oma=c(0,0,0,7))
  plot(x$metric[,1],main = main,ylim=ylim,xlab = xlab,ylab = ylab,type='l',col=cols[1])
  for (i in 2:9){
    lines(x$metric[,i],col=cols[i])
  }
  legend(length(x$vars)+length(x$vars)/20,0.5,yjust = 0.5,bty='n',lty=1,lwd=3,col = cols,
         legend = assesmentMetrics,xpd=NA)
}

plotCVtest <- function(x,main="",nam=NA,metrics=c("Accuracy","Kappa","F1","FDR","FNR","FPR")){
  if (any(is.na(nam))){
    nam <- names(x)
  }
  par(mfrow=c(2,3),oma=c(2,2,3,2))
  for (i in metrics){
    t <- sapply(x,function(y){y[i,]})
    boxplot(t,main=i,names = nam)
  }
  title(main,outer = TRUE)
}

## Plot ROC curve for cv result
rocCurve <- function(x,xlab='FPR',ylab='TPR',xlim=c(0,1),ylim=c(0,1),...){
  plot(x["FPR",],x["Sensitivity",],xlim = xlim,ylim=ylim,xlab = xlab,ylab = ylab,...)
  abline(0,1,lty=2)
}
