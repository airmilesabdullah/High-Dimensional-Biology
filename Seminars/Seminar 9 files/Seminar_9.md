Seminar 8
================
Abdullah Farouk
2018-03-14

``` r
library(MASS)
library(dplyr)
library(magrittr)
library(reshape)
library(car)
library(limma)
library(e1071)
library(glmnet)
library(ROCR)
library(CMA)
library(lattice)
library(class)
library(RCurl)
options('download.file.method'='curl')
library(GEOquery)
load("class_LNstatus.Rdata")
```

Data
====

``` r
# datgeo <- getGEO('GSE23177', GSEMatrix = TRUE) 
#   dat <- datgeo[[1]]   #Note that dat is an ExpressionSets
#   
#   str(pData(dat), max.level = 0)
#   
#   # extract only those variables of interest 
#   pData(dat) <-
#     subset(pData(dat),
#            select = c("characteristics_ch1.2",
#                       "characteristics_ch1.3","characteristics_ch1"))
#   names(pData(dat))<-c("LnStatus", "LnRatio", "Set")
# 
#   #Note: LNRatio will not be used in this Seminar. However, you can use it to try some of the regularization techniques learned in class
#   
#   # split the ExpressionSet into training and test sets. 
#   train.es <- dat[, dat$Set == "patient type: training set"]
#   test.es <- dat[ , dat$Set != "patient type: training set"]
# 
#   #Re-label factor
#   pData(train.es)$LnStatus <-
#       recode(pData(train.es)$LnStatus, "levels(pData(train.es)$LnStatus)[1]='neg'; else='pos'", levels = c('neg', 'pos'))
# 
#   pData(test.es)$LnStatus <-
#       recode(pData(test.es)$LnStatus, "levels(pData(test.es)$LnStatus)[1]='neg'; else='pos'",
#              levels = c('neg', 'pos'))
# 
#   # create data matrices with expression values (probesets in rows). Some of the functions we will use do not take ExpressionSets as objects
#   trainDat <- exprs(train.es)
#   testDat <- exprs(test.es)
# 
#   # Redefine the quantitative variable LnRatio to make it a numeric variable.
#   ntrain <- dim(pData(train.es))[1]
#   ntest <- dim(pData(test.es))[1]
#   
#   pData(train.es)$LnRatio <- as.numeric(unlist(strsplit(as.vector(unlist(pData(train.es)$LnRatio)), ":", fixed = TRUE))[(1:ntrain)*2])
#   pData(test.es)$LnRatio <- as.numeric(unlist(strsplit(as.vector(unlist(pData(test.es)$LnRatio)), ":", fixed = TRUE))[(1:ntest)*2])
# 
#   # save the data to avoid future re-downloading
#   save(dat,trainDat,testDat,train.es,test.es, file = "class_LNstatus.Rdata")
# ```
# 
# #Data Preprocessing
# ```{r}
# sum(is.na(trainDat))
# sum(is.na(trainDat))
```

Exercise 1: perform 100 runs of this CV before selecting a model to test! Add at least one rule to select the list of models, e.g., use genes with a p-val threshold &lt; cutoff.
=================================================================================================================================================================================

Here I changed the cutoff threshold to 0.9 Thus only the top 50 genes with p values lower than 0.9 are selected in each iteration.

``` r
set.seed(12)

reps <- 100
methods <- 7

cv_2 <- matrix(nrow = reps, ncol = methods)

for(j in 1:nrow(cv_2)) {
  nfold <- 6

tabTrain <- table(train.es$LnStatus)

indlist <- sapply(names(tabTrain), function(z) which(train.es$LnStatus == z), simplify = FALSE)

#Each row contains 8 pos and 8 negative samples. 

fold.pos <- matrix(sample(indlist[["pos"]]),nrow=nfold)
fold.neg <- matrix(sample(indlist[["neg"]]),nrow=nfold)

#Define here the constants that you will not evaluate. For example, I will use the top-50 limma genes

ngenes <- 50
nmethod <- 7 #number of methods you plan to compare. 

#Define here an output objects to store results
pr.err <- matrix(-1, nfold,nmethod, dimnames=list(paste0("Fold",1:nfold),c("1NN","5NN","10NN", "15NN","LDA","Logit","SVM")))

for(i in 1:nfold){

  #Test Fold for the i-th step
  testdat.fold<-trainDat[,c(fold.pos[i,],fold.neg[i,])]
  #I will create a factor of classes for the test set of the i_th fold
  testclass.fold<-train.es$LnStatus[c(fold.pos[i,],fold.neg[i,])]
  
    
  #The rest of the samples are the training set for the i-th step
  traindat.fold<-trainDat[,-c(fold.pos[i,],fold.neg[i,])]
  trainclass.fold<-train.es$LnStatus[-c(fold.pos[i,],fold.neg[i,])]

  #Step 1: feature selection (do you remember limma?). 

  # Note that a different set of genes will be selected for each fold! you can then compare how consistent these sets were.

  limma.dat<-as.data.frame(traindat.fold)
  desMat <- model.matrix(~ trainclass.fold, limma.dat) #design matrix
  trainFit <- lmFit(limma.dat, desMat)
  eBtrainFit <- eBayes(trainFit)
  
  # top-50 limma genes with a pvalue greater than 0.9
  top.fold <- topTable(eBtrainFit, coef = which(colnames(coef(trainFit)) != "(Intercept)"),
                       n = ngenes,sort.by="P", p.value = 0.9)
  
  #Retain the top-50 limma genes from the train and test sets
  traindat.fold <- traindat.fold[rownames(top.fold),]
  testdat.fold <-  testdat.fold[rownames(top.fold),]

  
  #STEP 2: select a classifier
  #Set a counter for the method tested
  l <- 0

  #kNN classifiers
  for(kk in c(1,5,10,15)) {
    #every time you get inside this loop, the l counter gets redefined (i.e., 1, 2, etc for         method 1, method 2, etc)
    l <- l+1

    #knn needs samples in rows
    yhat.knn <- knn(train=t(traindat.fold), test=t(testdat.fold), cl=trainclass.fold,
                    k = kk)
    #Store the prediction error for each kk within this fold
    pr.err[i,l]<- mean(testclass.fold != yhat.knn)
                          } #end of kNN loop

  #LDA method. Note that you can change the prior parameter to reflect a different proportion of case and control samples. The default is to use the class proportions from the training set.
  
  m.lda <- lda(x=t(traindat.fold), group=trainclass.fold, prior=c(.5, .5))
  yhat.lda <- predict(m.lda, newdata=t(testdat.fold))$class
  pr.err[i,"LDA"] <-mean(testclass.fold != yhat.lda)
   
  #Logit
  glm.dat <- data.frame(t(traindat.fold), group=trainclass.fold)
  
  # 50 factors still will cause optimization warnings  
  # Try without warning suppression to see 
  # To further reduce parameters, regularized regression can be used 
  # To use regularized regression uncomment lines followed by "uncomment for regularized regression" 
  suppressWarnings( m.log <- glm(group ~ ., data=glm.dat,family=binomial) ) 
  
  # uncomment for regularized regression 
  # m.log <- glmnet( t(traindat.fold) , trainclass.fold ,family="binomial") 

  pr.log <- predict(m.log,newdata=data.frame(t(testdat.fold)),type="response")
  
  # uncomment for regularized regression 
  # pr.log <- predict(m.log,newdata=data.frame(t(testdat.fold)),type="response",newx=t(testdat.fold)) 
  
  pr.cl <- rep(0,length(testclass.fold))
  pr.cl[pr.log > 1/2] <- "pos"
  pr.cl[pr.log <= 1/2] <- "neg"

  pr.cl <- factor(pr.cl)
  pr.err[i,"Logit"] <- mean( pr.cl != testclass.fold )

  #SVM
  m.svm <- svm(x=t(traindat.fold), y=trainclass.fold, cost=1, type="C-classification", 
               kernel="linear")
  pr.svm <- predict(m.svm,newdata=t(testdat.fold)) 
   
  pr.err[i,"SVM"] <- mean( pr.svm != testclass.fold )
  } #end of CV loop

  cv.err <- colMeans(pr.err)

  cv_2[j, ] <- rbind(cv.err)
  
}
colnames(cv_2) <- c("1NN", "5NN", "10NN", "15NN", "LDA", 'Logit', "SVM")
colMeans(cv_2)
```

    ##       1NN       5NN      10NN      15NN       LDA     Logit       SVM 
    ## 0.4082292 0.3843750 0.3847917 0.3788542 0.4607292 0.4630208 0.4302083
