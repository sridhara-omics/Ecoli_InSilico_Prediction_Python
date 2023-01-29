---
title: "FBA-GLMNET"
author: "Vish"
date: "2023-01-28"
output:
  pdf_document: default
  word_document: default
  html_document: default
---



Excerpt taken directly from publication: We initially carried out simulations on 49 growth conditions consisting of all pairwise combinations of 7 carbon and 7 nitrogen sources (Table 1). The compounds were selected from among the 174 carbon and 78 nitrogen sources previously used in Feist et al. The sources were chosen to yield distinct flux profiles, as assessed by k-means clustering of steady-state fluxes obtained using all pair-wise combinations of the 174 carbon and 78 nitrogen sources. Of the 7 carbon sources chosen, none resulted in any growth when used as a growth condition in the absence of a nitrogen source. By contrast, all nitrogen sources except ammonia yielded growth when supplied in the absence of a carbon source. Note that all of these remaining nitrogen sources do contain carbon atoms, thus additional supply of carbon was not strictly necessary for growth on these substrates.



```r
#' Set very small numbers to zero

threshold <- function(x, eps) {
  x[abs(x) < eps]=0
  return(x)
}
```


```r
#' Logical factor return to split data into train and test sets

split_train_test <- function(x, ratio){
  ## set seed for reproducibility
  set.seed(123)

  ## dplyr function
  sample <- sample(c(TRUE, FALSE), nrow(x), replace=TRUE, prob=c(1 - ratio, ratio))

  return(sample)

}
```


```r
#' K-fold cross validation of glmnet

runGLMNET <- function(inputs, targets, family, alpha, nfolds) {
  inputs = as.matrix(inputs)
  targets = as.factor(targets)

  cvGLMnet<-cv.glmnet(inputs, targets, family=family, alpha=alpha, standardize=TRUE, maxit=900000, type.measure="class", nfolds=nfolds, lambda.min.ratio=0.0002, nlambda=2000)
  return(cvGLMnet)
}
```


```r
#' Predict growth sources using trained model and test data

predictGLMNET <- function(trained_model, input_test){

  # derive coefficients from model
  s <- trained_model$lambda.min
  glmnet_coef <-as.vector(coef(trained_model, s=s))

  # predict using test data
  predicted_sources <- predict(trained_model, newx=as.matrix(input_test), type='class', s=s)

  return(predicted_sources)

}
```


```r
#' Confusion matrix - y versus yhat

confusion_matrix <- function(predicted_sources, actual_sources, number_sources){

  # fill out zeros for missing sources in the actual and predicted list
  pred_levels<-factor(as.numeric(paste(predicted_sources)),levels=c(1:number_sources))
  actual_levels<-factor(as.numeric(paste(actual_sources)),levels=c(1:number_sources))

  # confusion matrix to show actual values as columns and predicted as rows
  table_matrix<-table("Predicted"=pred_levels,"Actual"= actual_levels)

  return(table_matrix)

}
```


```r
#' Misclassification rate - yhat not predicted as y

misclassification_rate <- function(confusion_table) {

  rate <- 100 * (sum(confusion_table)-sum(diag(confusion_table))) / sum(confusion_table)

  return(rate)

}
```


```r
library(tidyverse)
library(glmnet)

# Read Flux Balance Analysis output file
M<-read.csv("C:/Users/Sridhara/Downloads/Ecoli_FBA_input_prediction-master/Ecoli_FBA_input_prediction-master/Analysis/RawData/FluxData49ReplicatesNoiseLevel1.csv",header=F)

dim(M)
```

```
## [1] 4900 2385
```
The FBA file has 4900 rows, corresponding to 100 replicates of 49 combinations of growth conditions. Out of 2385 rows, the first 2382 rows correspond to the simulated flux data, while the last 3 encode the target variable i.e, Carbon growth # (1-7), Nitrogen growth # (1-7) and Carbon-Nitrogen combination # (1-49).


```r
M[23:25, 2380:2385]
```

```
##        V2380 V2381     V2382 V2383 V2384 V2385
## 23 0.0015232     0 0.0015232     4     2    23
## 24 0.0044139     0 0.0044139     4     3    24
## 25 0.0018243     0 0.0018243     4     4    25
```
In the above table, an excerpt from table for columns through 2380:2385 shows simulated fluxes in the first 3 columns, while the last 3 columns shows the actual growth sources. We will use the simulated fluxes to train a model and then predict the actual growth conditions. In this markdown, we will predict only the combination i.e., V2385 column.



```r
# split the data into half for train and test
split_index <- split_train_test(M, 0.5)

Mtrain <- M[split_index, ]
Mtest <- M[!split_index, ]

col_index <- ncol(Mtrain)-3
inputs_train <- Mtrain[,1:col_index]
inputs_test <- Mtest[,1:col_index]

dim(inputs_train)
```

```
## [1] 2430 2382
```


```r
targets_train <- Mtrain[,ncol(Mtrain)]
targets_test <- Mtest[,ncol(Mtest)]
targets_train[1:5]
```

```
## [1] 2 4 5 7 8
```


```r
# Use alpha=1.0 for LASSO mode, nfolds = 3 (cv parameter) - Below is a k-crossfold validation of glmnet
# type.measure="class" applies to binomial and multinomial logistic regression only
glmnet_output <- runGLMNET(inputs_train, targets_train, "multinomial", 1.0, 3)
```

Text taken from cv.glmnet tutorial:
“class” gives misclassification error.
The multinomial model extends the binomial when the number of classes is more than two.
cv.glmnet has its special parameters including nfolds (the number of folds), foldid (usersupplied folds), and type.measure(the loss used for cross-validation).



```r
plot(glmnet_output)
```

![](Glmnet_testing_files/figure-latex/misclassificationError.pdf)<!-- --> 
Above plot is the cross-validation curve (red dotted line) along with upper and lower standard deviation curves
along the λ sequence (error bars).

Two special values along the λ sequence are indicated by the vertical
dotted lines. lambda.min is the value of λ that gives minimum mean cross-validated error, while lambda.1se
is the value of λ that gives the most regularized model such that the cross-validated error is within one
standard error of the minimum.



```r
# derive coefficients from the model to find the key features to predict the targets
# parameter s can be used to find the coefficients for a particular lambda
lasso_coefficients<-as.vector(coef(glmnet_output, s = glmnet_output$lambda.min))
```

The user can use the coefficients to find the key reactions that are critical to predict a particular growth source.


```r
# predict using test data
predicted_sources <- predictGLMNET(glmnet_output, inputs_test)
```

The metrics on the prediction accuracy can be calculated using the confusion matrix (actual versus predicted table) and/or the misclassification rate.


```r
# actual growth sources, and total number of combinations
actual_sources <- Mtest[,ncol(Mtest)]
number_sources <- length(unique(M[,ncol(M)]))

# confusion table
actual_vs_predicted <- confusion_matrix(predicted_sources, actual_sources, number_sources)
actual_vs_predicted
```

```
##          Actual
## Predicted  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23
##        1  49  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        2   0 46  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        3   0  0 52  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        4   0  0  0 44  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        5   0  0  0  0 59  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        6   0  0  0  0  0 46  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        7   0  0  0  0  0  0 48  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        8   0  0  0  0  0  0  0 50  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0
##        9   0  0  0  0  0  0  0  0 46  0  0  0  0  0  0  0  0  0  0  0  0  0  1
##        10  0  0  0  0  0  0  0  0  0 45  0  0  0  0  0  0  0  0  0  0  0  0  1
##        11  0  0  0  0  0  0  0  0  0  0 55  0  0  0  0  0  0  0  0  0  0  0  0
##        12  0  0  0  0  0  0  0  0  0  0  0 54  0  0  0  0  0  0  0  0  0  0  0
##        13  0  0  0  0  0  0  0  0  1  0  0  0 50  0  0  0  0  0  0  0  0  0  0
##        14  0  0  0  0  0  0  0  0  0  0  0  0  0 43  0  0  0  0  0  0  0  0  0
##        15  0  0  0  0  0  0  0  0  0  0  0  0  0  0 64  2  0  0  0  0  1  0  0
##        16  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 52  0  0  0  0  0  0  0
##        17  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 50  0  0  0  0  0  0
##        18  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 50  0  0  0  0  0
##        19  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 45  0  0  0  0
##        20  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0 46  0  0  0
##          Actual
## Predicted 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46
##        1   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        2   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        3   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        4   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        5   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        6   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        7   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0
##        8   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        9   0  0  0  0  0  0  1  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0
##        10  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        11  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        12  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
##        13  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        14  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        15  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        16  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        17  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        18  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        19  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##        20  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
##          Actual
## Predicted 47 48 49
##        1   0  0  0
##        2   0  0  1
##        3   0  0  0
##        4   0  0  0
##        5   0  0  0
##        6   0  0  0
##        7   0  0  0
##        8   0  0  0
##        9   0  0  0
##        10  0  0  0
##        11  0  0  0
##        12  0  0  0
##        13  0  0  0
##        14  0  0  0
##        15  0  0  0
##        16  0  0  0
##        17  0  0  0
##        18  0  0  0
##        19  0  0  0
##        20  0  0  0
##  [ reached getOption("max.print") -- omitted 29 rows ]
```
The diagonal of the confusion table shows the correctly predicted sources (i.e., 1 vs 1, 2 vs 2, 3 vs 3). the off-diagonal numbers show the misclassified sources. For example, one each of 12th and 14th combinations of the growth source is wrongly predicted as 8th source. However the correctly classified numbers for 12th and 14th are 53 and 43 respectively.

Another way to present the misclassified data is by calculating the misclassification rate using the confusion table.

```r
# misclassification rate
misclassification_rate_CN <- misclassification_rate(actual_vs_predicted)
misclassification_rate_CN
```

```
## [1] 2.42915
```
The misclassification rate is 2.388%, which confirms that the model has a very good accuracy of predicting the growth conditions of Carbon and Nitrogen combinations. The same analyses can be extended to do different predictions (i.e, separate predictions of carbon (C) and nitrogen (N)).
