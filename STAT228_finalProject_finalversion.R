##########################################
#
# STAT 228 Final Project Code
# Version Date: May 2, 2023
#
##########################################

########################################## Libraries ##########################################
library("usdm")
library("glmulti")
library("nnet")
library("tree") # CART
library("e1071") # SVM
library("randomForest") # bagging
library("lmtest") # Maximum Likelihood Test

########################################## Data ##########################################
# codon.data <- read.csv("codon_usage.csv", header = TRUE)
codon.data <- read.csv("/Users/alicezhang/Desktop/Spring 2023/stat228_rdata/codon_usage.csv", header = TRUE)
table(unique(codon.data$Kingdom))
table((codon.data$Kingdom))
dim(codon.data)

### data cleaning
# Clean the Kingdom column by removing leading and trailing whitespace
codon.data$Kingdom <- trimws(codon.data$Kingdom)
# Convert all values in the Kingdom column to lowercase
codon.data$Kingdom <- tolower(codon.data$Kingdom)
# Check the unique values in the Kingdom column to ensure that there are no formatting differences
unique(codon.data$Kingdom)

### subset data, mam	pri	rod	vrt
data.sub <- subset(codon.data, codon.data$Kingdom %in% c("mam", "vrt", "pri", "rod"))
data.sub <- data.sub[,-c(2,3,4,5)]
dim(data.sub)
data.sub[, 2:65] <- apply(data.sub[, 2:65], 2, as.numeric)
head(data.sub)
data.sub$Kingdom <- as.factor(data.sub$Kingdom)
X <- data.sub[, 2:65]
levels(data.sub$Kingdom)

############################## multicollin w/ vifstep ########################################

vifstep(X, th=10)
X = X[, !colnames(X) %in% c("CUA", "UGA", "GAG", "AAG", "CAG")]
data.sub.nomc = data.sub[, !colnames(data.sub) %in% c("CUA", "UGA", "GAG", "AAG", "CAG")]
colnames(data.sub.nomc)
dim(data.sub.nomc)

###############################################################################
## Cross Validation Preparation
##
## Partition Data 8:2 (training, validation)
###############################################################################

set.seed(228)
n <- dim(data.sub.nomc)[1]
n.shuffle <- sample(1:n, n, replace = FALSE) # shuffle the n indices
data.train = data.sub.nomc[n.shuffle[1:2435],]	# 80% of 3044 is roughly 2435
data.validate = data.sub.nomc[n.shuffle[2436:n],]

###############################################################################
## Variable Selection + Model Building
##
## Multinomial Ordinal Regression
##
## Stepwise Selection based on AIC and BIC
###############################################################################

fit.multinomial.full <- multinom(Kingdom~., data=data.train)
summary(fit.multinomial.full)

aic.model <- step(fit.multinomial.full, direction="both", trace=0,k=2) # final aic value 934.922509 
summary(aic.model)
bic.model <- step(fit.multinomial.full, direction="both", trace=0,k=log(n)) # final bic value
summary(bic.model)

# AIC model
aic.data.train <- data.train[,-c(16, 22, 30, 31, 39, 41, 52, 55, 60)]
aic.data.validate <- data.validate[,-c(16, 22, 30, 31, 39, 41, 52, 55, 60)]
fit.multinomial.AIC <- multinom(Kingdom~., data=aic.data.train)

# BIC model
BIC.data.train <- data.train[,c(1,2,3,4,6,9,11:15, 17:21,23,27,28,32, 37,40, 42,44,45,49,50,51,56, 58,59)]
BIC.data.validate <- data.validate[,c(1,2,3,4,6,9,11:15, 17:21,23,27,28,32, 37,40, 42,44,45,49,50,51,56, 58,59)]
fit.multinomial.BIC <- multinom(Kingdom~., data=BIC.data.train)

# Perform the LRT
lrt <- lrtest(fit.multinomial.AIC, fit.multinomial.BIC)

########## Model 1: Stepwise AIC ########## 

## 1. Performance Metrics (accuracy, sensitivity, specificity, precision)
# Partition the data into k equal subsets
set.seed(228)
K <- 10 # 10-fold CV
n.fold <- floor((dim(data.validate)[1])/K)
n.shuffle <- sample(1:(dim(data.validate)[1]), (dim(data.validate)[1]), replace = FALSE) # shuffle the n indices
index.fold <- list() 

for (i in 1:K){
  if (i < K){
    index.fold[[i]] <- n.shuffle[((i-1) * n.fold + 1):(i * n.fold)]
  } 
  else{
    index.fold[[i]] <- n.shuffle[((K-1) * n.fold + 1):(dim(data.validate)[1])]
  }
}

CV_misclassification.AIC <- NA
CV_sensitivity_AIC <- NA
CV_specificity_AIC <- NA
CV_precision_AIC <- NA


for (kval in 1:K) {
  fit <- multinom(Kingdom ~ .,
                  data = aic.data.validate[-index.fold[[kval]], ])
  
  Y.hat <-
    predict(fit,  newdata = aic.data.validate[index.fold[[kval]], ], type = "class") #<- ifelse(p.hat > c[i], 1, 0)
  confusion.table <-
    table(Y.hat, aic.data.validate[index.fold[[kval]], ]$Kingdom)
  
  CV_misclassification.AIC <-
    (
      confusion.table[1, 2] + confusion.table[1, 3] + confusion.table[1, 4]
      + confusion.table[2, 1] + confusion.table[2, 3] +
        confusion.table[2, 4]
      + confusion.table[3, 1] + confusion.table[3, 2] +
        confusion.table[2, 4]
      + confusion.table[4, 1] + confusion.table[4, 2] +
        confusion.table[4, 3]
    ) / dim(BIC.data.validate[index.fold[[kval]],])[1]
  
  ### sensitivity
  sensitivity.AIC <- rep(0, 4)
  
  for (ii in 1:4) {
    # print(i)
    if (sum(confusion.table[, ii]) == 0) {
      sensitivity.AIC[i] <- 0.0
    }
    else{
      sensitivity.AIC[ii] <-
        confusion.table[, ii][ii] / sum(confusion.table[, ii])
    }
  }
  
  CV_sensitivity_AIC <-
    (
      sensitivity.AIC[1] * sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "mam") + sensitivity.AIC[2] *
        sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "pri") + sensitivity.AIC[3] *
        sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "rod") + sensitivity.AIC[4] *
        sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "vrt")
    ) / dim(aic.data.validate[index.fold[[kval]],])[1]
  
  ### specificity = true neg / (true neg + false pos)
  specificity.AIC <- rep(0, 4)
  TN.mam <-
    confusion.table[2, 2] + confusion.table[2, 3] + confusion.table[2, 4] + confusion.table[3, 2] +
    confusion.table[3, 3] + confusion.table[3, 4] + confusion.table[4, 2] + confusion.table[4, 3] +
    confusion.table[4, 4]
  FP.mam <-
    confusion.table[1, 2] + confusion.table[1, 3] + confusion.table[1, 4]
  specificity.mam <- TN.mam / (TN.mam + FP.mam)
  specificity.AIC[1] <- specificity.mam
  
  TN.pri <-
    confusion.table[1, 1] + confusion.table[1, 3] + confusion.table[1, 4] + confusion.table[3, 1] +
    confusion.table[3, 3] + confusion.table[3, 4] + confusion.table[4, 1] + confusion.table[4, 3] +
    confusion.table[4, 4]
  FP.pri <-
    confusion.table[2, 1] + confusion.table[2, 3] + confusion.table[2, 4]
  specificity.pri <- TN.pri / (TN.pri + FP.pri)
  specificity.AIC[2] <- specificity.pri
  
  TN.rod <-
    confusion.table[1, 1] + confusion.table[1, 2] + confusion.table[1, 4] + confusion.table[2, 1] +
    confusion.table[2, 2] + confusion.table[2, 4] + confusion.table[4, 1] + confusion.table[4, 2] +
    confusion.table[4, 4]
  FP.rod <-
    confusion.table[3, 1] + confusion.table[3, 2] + confusion.table[3, 4]
  specificity.rod <- TN.rod / (TN.rod + FP.rod)
  specificity.AIC[3] <- specificity.rod
  
  TN.vrt <-
    confusion.table[1, 1] + confusion.table[1, 2] + confusion.table[1, 3] + confusion.table[2, 1] +
    confusion.table[2, 2] + confusion.table[2, 3] + confusion.table[3, 1] + confusion.table[3, 2] + confusion.table[3, 3]
  FP.vrt <-
    confusion.table[4, 1] + confusion.table[4, 2] + confusion.table[4, 3]
  specificity.vrt <- TN.vrt / (TN.vrt + FP.vrt)
  specificity.AIC[4] <- specificity.vrt
  
  CV_specificity_AIC <-
    (
      specificity.AIC[1] * sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "mam") + specificity.AIC[2] *
        sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "pri") + specificity.AIC[3] *
        sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "rod") + specificity.AIC[4] *
        sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "vrt")
    ) / dim(aic.data.validate[index.fold[[kval]],])[1]
  
  
  ### precision
  precision.AIC <- rep(0, 4)
  for (ii in 1:4) {
    # print(i)
    if (sum(confusion.table[ii, ]) == 0) {
      precision.AIC[ii] <- 0.0
    }
    else{
      precision.AIC[ii] <-
        confusion.table[ii, ][ii] / sum(confusion.table[ii, ])
    }
  }
  CV_precision_AIC <-
    (
      precision.AIC[1] * sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "mam") + precision.AIC[2] *
        sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "pri") + precision.AIC[3] *
        sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "rod") + precision.AIC[4] *
        sum(aic.data.validate[index.fold[[kval]],]$Kingdom == "vrt")
    ) / dim(aic.data.validate[index.fold[[kval]],])[1]
  
  
}




fmeasure.AIC <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (precision.AIC[i]+sensitivity.AIC[i] == 0){
    fmeasure.AIC[i] <- 0.0
  }
  else{
    fmeasure.AIC[i] <- 2*(precision.AIC[i]*sensitivity.AIC[i])/(precision.AIC[i]+sensitivity.AIC[i])
  }
}

### Performance Metrics
sensitivity.AIC
specificity.AIC
precision.AIC
fmeasure.AIC

### Aggregate Performance Metrics
perf.metrics_AIC_wa <- data.frame(cbind(CV_misclassification.AIC, CV_sensitivity_AIC, 
                                     CV_specificity_AIC, CV_precision_AIC))
perf.metrics_AIC_wa$F.measure <- 2 * CV_precision_AIC * CV_sensitivity_AIC / (CV_precision_AIC + CV_sensitivity_AIC)
perf.metrics_AIC_wa

########## Model 2: Stepwise BIC ########## 

## 1. Performance Metrics (accuracy, sensitivity, specificity, precision)

CV_misclassification.BIC <- NA
CV_sensitivity_BIC <- NA
CV_specificity_BIC <- NA
CV_precision_BIC <- NA


for (kval in 1:K) {
  fit <- multinom(Kingdom ~ .,
                  data = BIC.data.validate[-index.fold[[kval]], ])
    Y.hat <-
    predict(fit,  newdata = BIC.data.validate[index.fold[[kval]], ], type = "class") #<- ifelse(p.hat > c[i], 1, 0)
  confusion.table <-
    table(Y.hat, BIC.data.validate[index.fold[[kval]], ]$Kingdom)
  
  CV_misclassification.BIC <-
    (
      confusion.table[1, 2] + confusion.table[1, 3] + confusion.table[1, 4]
      + confusion.table[2, 1] + confusion.table[2, 3] +
        confusion.table[2, 4]
      + confusion.table[3, 1] + confusion.table[3, 2] +
        confusion.table[2, 4]
      + confusion.table[4, 1] + confusion.table[4, 2] +
        confusion.table[4, 3]
    ) / dim(BIC.data.validate[index.fold[[kval]],])[1]
  
  ### Sensitivity
  sensitivity.BIC <- rep(0, 4)
  
  for (ii in 1:4) {
    # print(i)
    if (sum(confusion.table[, ii]) == 0) {
      sensitivity.BIC[i] <- 0.0
    }
    else{
      sensitivity.BIC[ii] <-
        confusion.table[, ii][ii] / sum(confusion.table[, ii])
    }
  }
  
  CV_sensitivity_BIC <-
    (
      sensitivity.BIC[1] * sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "mam") + sensitivity.BIC[2] *
        sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "pri") + sensitivity.BIC[3] *
        sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "rod") + sensitivity.BIC[4] *
        sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "vrt")
    ) / dim(BIC.data.validate[index.fold[[kval]],])[1]
  
  ### Specificity = true neg / (true neg + false pos)
  specificity.BIC <- rep(0, 4)
  TN.mam <-
    confusion.table[2, 2] + confusion.table[2, 3] + confusion.table[2, 4] + confusion.table[3, 2] +
    confusion.table[3, 3] + confusion.table[3, 4] + confusion.table[4, 2] + confusion.table[4, 3] +
    confusion.table[4, 4]
  FP.mam <-
    confusion.table[1, 2] + confusion.table[1, 3] + confusion.table[1, 4]
  specificity.mam <- TN.mam / (TN.mam + FP.mam)
  specificity.BIC[1] <- specificity.mam
  
  TN.pri <-
    confusion.table[1, 1] + confusion.table[1, 3] + confusion.table[1, 4] + confusion.table[3, 1] +
    confusion.table[3, 3] + confusion.table[3, 4] + confusion.table[4, 1] + confusion.table[4, 3] +
    confusion.table[4, 4]
  FP.pri <-
    confusion.table[2, 1] + confusion.table[2, 3] + confusion.table[2, 4]
  specificity.pri <- TN.pri / (TN.pri + FP.pri)
  specificity.BIC[2] <- specificity.pri
  
  TN.rod <-
    confusion.table[1, 1] + confusion.table[1, 2] + confusion.table[1, 4] + confusion.table[2, 1] +
    confusion.table[2, 2] + confusion.table[2, 4] + confusion.table[4, 1] + confusion.table[4, 2] +
    confusion.table[4, 4]
  FP.rod <-
    confusion.table[3, 1] + confusion.table[3, 2] + confusion.table[3, 4]
  specificity.rod <- TN.rod / (TN.rod + FP.rod)
  specificity.BIC[3] <- specificity.rod
  
  TN.vrt <-
    confusion.table[1, 1] + confusion.table[1, 2] + confusion.table[1, 3] + confusion.table[2, 1] +
    confusion.table[2, 2] + confusion.table[2, 3] + confusion.table[3, 1] + confusion.table[3, 2] + confusion.table[3, 3]
  FP.vrt <-
    confusion.table[4, 1] + confusion.table[4, 2] + confusion.table[4, 3]
  specificity.vrt <- TN.vrt / (TN.vrt + FP.vrt)
  specificity.BIC[4] <- specificity.vrt
  
  CV_specificity_BIC <-
    (
      specificity.BIC[1] * sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "mam") + specificity.BIC[2] *
        sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "pri") + specificity.BIC[3] *
        sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "rod") + specificity.BIC[4] *
        sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "vrt")
    ) / dim(BIC.data.validate[index.fold[[kval]],])[1]
  
  
  ### Precision
  precision.BIC <- rep(0, 4)
  for (ii in 1:4) {
    # print(i)
    if (sum(confusion.table[ii, ]) == 0) {
      precision.BIC[ii] <- 0.0
    }
    else{
      precision.BIC[ii] <-
        confusion.table[ii, ][ii] / sum(confusion.table[ii, ])
    }
  }
  CV_precision_BIC <-
    (
      precision.BIC[1] * sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "mam") + precision.BIC[2] *
        sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "pri") + precision.BIC[3] *
        sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "rod") + precision.BIC[4] *
        sum(BIC.data.validate[index.fold[[kval]],]$Kingdom == "vrt")
    ) / dim(BIC.data.validate[index.fold[[kval]],])[1]
  
  
}

fmeasure.BIC <- rep(0, 4)
for (i in 1:4){
  # print(i)
  if (precision.BIC[i]+sensitivity.BIC[i] == 0){
    fmeasure.BIC[i] <- 0.0
  }
  else{
    fmeasure.BIC[i] <- 2*(precision.BIC[i]*sensitivity.BIC[i])/(precision.BIC[i]+sensitivity.BIC[i])
  }
}

### Performance Metrics
sensitivity.BIC
specificity.BIC
precision.BIC
fmeasure.BIC

### Aggregate Performance Metrics
perf.metrics_BIC_wa <- data.frame(cbind(CV_misclassification.AIC, CV_sensitivity_AIC, 
                                     CV_specificity_AIC, CV_precision_AIC))
perf.metrics_BIC_wa$F.measure <- 2 * CV_precision_AIC * CV_sensitivity_AIC / (CV_precision_AIC + CV_sensitivity_AIC)
perf.metrics_BIC_wa

###############################################################################
## 
## CART, RF, Bagging
## 
###############################################################################

########## Model 3:CART ########## 

### 1. Performance Metrics (accuracy, sensitivity, specificity, precision)

set.seed(228)
tree.cart <- tree(as.factor(Kingdom)~., data=data.train)
plot(tree.cart)
text(tree.cart, pretty=0)
result.cart <- cv.tree(tree.cart, K=10, FUN=prune.tree)
plot(result.cart)
# optimal size = 14
tree.new.cart = prune.tree(tree.cart, best=14)
plot(tree.new.cart)
text(tree.new.cart, pretty=0)
tree.new.cart
summary(tree.new.cart)

Y.cart.hat <- predict(tree.new.cart, newdata=data.validate, type="class")
length(Y.cart.hat)

### misclassification
confusion.cart <- table(Y.cart.hat, data.validate$Kingdom)
misclassification.cart <- (confusion.cart[1,2]+confusion.cart[1,3]+confusion.cart[1,4]
                           +confusion.cart[2,1]+confusion.cart[2,3]+confusion.cart[2,4]
                           +confusion.cart[3,1]+confusion.cart[3,2]+confusion.cart[2,4]
                           +confusion.cart[4,1]+confusion.cart[4,2]+confusion.cart[4,3])/length(Y.cart.hat)

### sensitivity
sensitivity.cart <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (sum(confusion.cart[,i]) == 0){
    sensitivity.cart[i] <- 0.0
  }
  else{
    sensitivity.cart[i] <- confusion.cart[,i][i]/sum(confusion.cart[,i])
  }
}

sensitivity.cart.wa <-
  (
    sensitivity.cart[1] * sum(data.validate$Kingdom == "mam") + sensitivity.cart[2] *
      sum(data.validate$Kingdom == "pri") + sensitivity.cart[3] *
      sum(data.validate$Kingdom == "rod") + sensitivity.cart[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]

# 1st element in vector is for mam, 2nd is for pri, 3rd is for rod, 4th is for vrt

### specificity = true neg / (true neg + false pos)
TN.mam <- confusion.cart[2, 2] + confusion.cart[2, 3] + confusion.cart[2, 4] + confusion.cart[3, 2] +
  confusion.cart[3, 3] + confusion.cart[3, 4] + confusion.cart[4, 2] + confusion.cart[4, 3] +
  confusion.cart[4, 4]
FP.mam <- confusion.cart[1,2]+confusion.cart[1,3]+confusion.cart[1,4]
specificity.mam <- TN.mam / (TN.mam+FP.mam)

TN.pri <-
  confusion.cart[1, 1] + confusion.cart[1, 3] + confusion.cart[1, 4] + confusion.cart[3, 1] +
  confusion.cart[3, 3] + confusion.cart[3, 4] + confusion.cart[4, 1] + confusion.cart[4, 3] +
  confusion.cart[4, 4]
FP.pri <- confusion.cart[2, 1] + confusion.cart[2, 3] + confusion.cart[2, 4] 
specificity.pri <- TN.pri / (TN.pri+FP.pri)

TN.rod <- confusion.cart[1, 1] + confusion.cart[1, 2] + confusion.cart[1, 4] + confusion.cart[2, 1] +
  confusion.cart[2, 2] + confusion.cart[2, 4] + confusion.cart[4, 1] + confusion.cart[4, 2] +
  confusion.cart[4, 4]
FP.rod <- confusion.cart[3, 1] + confusion.cart[3, 2] + confusion.cart[3, 4]
specificity.rod <- TN.rod / (TN.rod+FP.rod)

TN.vrt <-
  confusion.cart[1, 1] + confusion.cart[1, 2] + confusion.cart[1, 3] + confusion.cart[2, 1] +
  confusion.cart[2, 2] + confusion.cart[2, 3] + confusion.cart[3, 1] + confusion.cart[3, 2] + confusion.cart[3, 3]
FP.vrt <- confusion.cart[4, 1] + confusion.cart[4, 2] + confusion.cart[4, 3]
specificity.vrt <- TN.vrt / (TN.vrt+FP.vrt)
specificity.cart <- c(specificity.mam, specificity.pri, specificity.rod, specificity.vrt)

specificity.cart.wa <-
  (
    specificity.cart[1] * sum(data.validate$Kingdom == "mam") + specificity.cart[2] *
      sum(data.validate$Kingdom == "pri") + specificity.cart[3] *
      sum(data.validate$Kingdom == "rod") + specificity.cart[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]

### precision
precision.cart <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (sum(confusion.cart[i,]) == 0){
    precision.cart[i] <- 0.0
  }
  else{
    precision.cart[i] <- confusion.cart[i,][i]/sum(confusion.cart[i,])
  }
}

precision.cart.wa <-
  (
    precision.cart[1] * sum(data.validate$Kingdom == "mam") + precision.cart[2] *
      sum(data.validate$Kingdom == "pri") + precision.cart[3] *
      sum(data.validate$Kingdom == "rod") + precision.cart[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]

### f-measure
fmeasure.cart <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (precision.cart[i]+sensitivity.cart[i] == 0){
    fmeasure.cart[i] <- 0.0
  }
  else{
    fmeasure.cart[i] <- 2*(precision.cart[i]*sensitivity.cart[i])/(precision.cart[i]+sensitivity.cart[i])
  }
}
fmeasure.cart


perf.metrics_cart <- data.frame(cbind(misclassification.cart, sensitivity.cart.wa, 
                                      specificity.cart.wa, precision.cart.wa))


### F-measure = (precision*recall) / (precision+recall)
perf.metrics_cart$F.measure <- 2 * precision.cart.wa * sensitivity.cart.wa / (precision.cart.wa + sensitivity.cart.wa)
perf.metrics_cart

########## Model 4:Bagging ########## 
## Bagging

tree.bagging <- randomForest(as.factor(Kingdom)~., data=data.train, mtry=dim(data.train)[2]-1)
tree.bagging
plot(tree.bagging)

Y.bagging.hat <- predict(tree.bagging, newdata=data.validate, type="class")
table(Y.bagging.hat, data.validate$Kingdom)
confusion.bagging <- table(Y.bagging.hat, data.validate$Kingdom)

### misclassification
misclassification.bagging <- (confusion.bagging[1,2]+ confusion.bagging[1,3]+ confusion.bagging[1,4]
                              +confusion.bagging[2,1]+ confusion.bagging[2,3]+ confusion.bagging[2,4]
                              + confusion.bagging[3,1]+ confusion.bagging[3,2]+ confusion.bagging[2,4]
                              + confusion.bagging[4,1]+ confusion.bagging[4,2]+ confusion.bagging[4,3])/length(Y.bagging.hat)
misclassification.bagging

sensitivity.bagging <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (sum(confusion.bagging[,i]) == 0){
    sensitivity.bagging[i] <- 0.0
  }
  else{
    sensitivity.bagging[i] <- confusion.bagging[,i][i]/sum(confusion.bagging[,i])
  }
}

sensitivity.bagging.wa <-
  (
    sensitivity.bagging[1] * sum(data.validate$Kingdom == "mam") + sensitivity.bagging[2] *
      sum(data.validate$Kingdom == "pri") + sensitivity.bagging[3] *
      sum(data.validate$Kingdom == "rod") + sensitivity.bagging[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]

# 1st element in vector is for mam, 2nd is for pri, 3rd is for rod, 4th is for vrt

### specificity = true neg / (true neg + false pos)
TN.mam <- confusion.bagging[2, 2] + confusion.bagging[2, 3] + confusion.bagging[2, 4] + confusion.bagging[3, 2] +
  confusion.bagging[3, 3] + confusion.bagging[3, 4] + confusion.bagging[4, 2] + confusion.bagging[4, 3] +
  confusion.bagging[4, 4]
FP.mam <- confusion.bagging[1,2]+confusion.bagging[1,3]+confusion.bagging[1,4]
specificity.mam <- TN.mam / (TN.mam+FP.mam)

TN.pri <-
  confusion.bagging[1, 1] + confusion.bagging[1, 3] + confusion.bagging[1, 4] + confusion.bagging[3, 1] +
  confusion.bagging[3, 3] + confusion.bagging[3, 4] + confusion.bagging[4, 1] + confusion.bagging[4, 3] +
  confusion.bagging[4, 4]
FP.pri <- confusion.bagging[2, 1] + confusion.bagging[2, 3] + confusion.bagging[2, 4] 
specificity.pri <- TN.pri / (TN.pri+FP.pri)

TN.rod <- confusion.bagging[1, 1] + confusion.bagging[1, 2] + confusion.bagging[1, 4] + confusion.bagging[2, 1] +
  confusion.bagging[2, 2] + confusion.bagging[2, 4] + confusion.bagging[4, 1] + confusion.bagging[4, 2] +
  confusion.bagging[4, 4]
FP.rod <- confusion.bagging[3, 1] + confusion.bagging[3, 2] + confusion.bagging[3, 4]
specificity.rod <- TN.rod / (TN.rod+FP.rod)

TN.vrt <-
  confusion.bagging[1, 1] + confusion.bagging[1, 2] + confusion.bagging[1, 3] + confusion.bagging[2, 1] +
  confusion.bagging[2, 2] + confusion.bagging[2, 3] + confusion.bagging[3, 1] + confusion.bagging[3, 2] + confusion.bagging[3, 3]
FP.vrt <- confusion.bagging[4, 1] + confusion.bagging[4, 2] + confusion.bagging[4, 3]
specificity.vrt <- TN.vrt / (TN.vrt+FP.vrt)
specificity.bagging <- c(specificity.mam, specificity.pri, specificity.rod, specificity.vrt)

specificity.bagging.wa <-
  (
    specificity.bagging[1] * sum(data.validate$Kingdom == "mam") + specificity.bagging[2] *
      sum(data.validate$Kingdom == "pri") + specificity.bagging[3] *
      sum(data.validate$Kingdom == "rod") + specificity.bagging[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]


precision.bagging <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (sum(confusion.bagging[i,]) == 0){
    precision.bagging[i] <- 0.0
  }
  else{
    precision.bagging[i] <- confusion.bagging[i,][i]/sum(confusion.bagging[i,])
  }
}

precision.bagging.wa <-
  (
    precision.bagging[1] * sum(data.validate$Kingdom == "mam") + precision.bagging[2] *
      sum(data.validate$Kingdom == "pri") + precision.bagging[3] *
      sum(data.validate$Kingdom == "rod") + precision.bagging[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]

fmeasure.bagging <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (precision.bagging[i]+sensitivity.bagging[i] == 0){
    fmeasure.bagging[i] <- 0.0
  }
  else{
    fmeasure.bagging[i] <- 2*(precision.bagging[i]*sensitivity.bagging[i])/(precision.bagging[i]+sensitivity.bagging[i])
  }
}
fmeasure.bagging


perf.metrics_bagging <- data.frame(cbind(misclassification.bagging, sensitivity.bagging.wa, 
                                         specificity.bagging.wa, precision.bagging.wa))


### F-measure = (precision*recall) / (precision+recall)
perf.metrics_bagging$F.measure <- 2 * precision.bagging.wa * sensitivity.bagging.wa / (precision.bagging.wa + sensitivity.bagging.wa)
perf.metrics_bagging

########## Model 5:Random Forest ########## 
## Random Forest
set.seed(228)
fit.RF <- randomForest(as.factor(Kingdom)~., data=data.train, mtry=floor(sqrt(dim(data.train)[2]-1)))
fit.RF
plot(fit.RF)
text(fit.RF, pretty=0)

Y.RF.hat <- predict(fit.RF, newdata=data.validate, type="class")
table(Y.RF.hat, data.validate$Kingdom)
confusion.RF <- table(Y.RF.hat, data.validate$Kingdom)

### misclassification
misclassification.RF <- (confusion.RF[1,2]+ confusion.RF[1,3]+ confusion.RF[1,4]
                         +confusion.RF[2,1]+ confusion.RF[2,3]+ confusion.RF[2,4]
                         + confusion.RF[3,1]+ confusion.RF[3,2]+ confusion.RF[2,4]
                         + confusion.RF[4,1]+ confusion.RF[4,2]+ confusion.RF[4,3])/length(Y.RF.hat)
misclassification.RF

sensitivity.RF <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (sum(confusion.RF[,i]) == 0){
    sensitivity.RF[i] <- 0.0
  }
  else{
    sensitivity.RF[i] <- confusion.RF[,i][i]/sum(confusion.RF[,i])
  }
}

sensitivity.RF.wa <-
  (
    sensitivity.RF[1] * sum(data.validate$Kingdom == "mam") + sensitivity.RF[2] *
      sum(data.validate$Kingdom == "pri") + sensitivity.RF[3] *
      sum(data.validate$Kingdom == "rod") + sensitivity.RF[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]

# 1st element in vector is for mam, 2nd is for pri, 3rd is for rod, 4th is for vrt

### specificity = true neg / (true neg + false pos)
TN.mam <- confusion.RF[2, 2] + confusion.RF[2, 3] + confusion.RF[2, 4] + confusion.RF[3, 2] +
  confusion.RF[3, 3] + confusion.RF[3, 4] + confusion.RF[4, 2] + confusion.RF[4, 3] +
  confusion.RF[4, 4]
FP.mam <- confusion.RF[1,2]+confusion.RF[1,3]+confusion.RF[1,4]
specificity.mam <- TN.mam / (TN.mam+FP.mam)

TN.pri <-
  confusion.RF[1, 1] + confusion.RF[1, 3] + confusion.RF[1, 4] + confusion.RF[3, 1] +
  confusion.RF[3, 3] + confusion.RF[3, 4] + confusion.RF[4, 1] + confusion.RF[4, 3] +
  confusion.RF[4, 4]
FP.pri <- confusion.RF[2, 1] + confusion.RF[2, 3] + confusion.RF[2, 4] 
specificity.pri <- TN.pri / (TN.pri+FP.pri)

TN.rod <- confusion.RF[1, 1] + confusion.RF[1, 2] + confusion.RF[1, 4] + confusion.RF[2, 1] +
  confusion.RF[2, 2] + confusion.RF[2, 4] + confusion.RF[4, 1] + confusion.RF[4, 2] +
  confusion.RF[4, 4]
FP.rod <- confusion.RF[3, 1] + confusion.RF[3, 2] + confusion.RF[3, 4]
specificity.rod <- TN.rod / (TN.rod+FP.rod)

TN.vrt <-
  confusion.RF[1, 1] + confusion.RF[1, 2] + confusion.RF[1, 3] + confusion.RF[2, 1] +
  confusion.RF[2, 2] + confusion.RF[2, 3] + confusion.RF[3, 1] + confusion.RF[3, 2] + confusion.RF[3, 3]
FP.vrt <- confusion.RF[4, 1] + confusion.RF[4, 2] + confusion.RF[4, 3]
specificity.vrt <- TN.vrt / (TN.vrt+FP.vrt)
specificity.RF <- c(specificity.mam, specificity.pri, specificity.rod, specificity.vrt)

specificity.RF.wa <-
  (
    specificity.RF[1] * sum(data.validate$Kingdom == "mam") + specificity.RF[2] *
      sum(data.validate$Kingdom == "pri") + specificity.RF[3] *
      sum(data.validate$Kingdom == "rod") + specificity.RF[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]


precision.RF <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (sum(confusion.RF[i,]) == 0){
    precision.RF[i] <- 0.0
  }
  else{
    precision.RF[i] <- confusion.RF[i,][i]/sum(confusion.RF[i,])
  }
}

precision.RF.wa <-
  (
    precision.RF[1] * sum(data.validate$Kingdom == "mam") + precision.RF[2] *
      sum(data.validate$Kingdom == "pri") + precision.RF[3] *
      sum(data.validate$Kingdom == "rod") + precision.RF[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]


## F Measure
fmeasure.RF <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (precision.RF[i]+sensitivity.RF[i] == 0){
    fmeasure.RF[i] <- 0.0
  }
  else{
    fmeasure.RF[i] <- 2*(precision.RF[i]*sensitivity.RF[i])/(precision.RF[i]+sensitivity.RF[i])
  }
}
fmeasure.RF


perf.metrics_RF <- data.frame(cbind(misclassification.RF, sensitivity.RF.wa, 
                                    specificity.RF.wa, precision.RF.wa))

### F-measure = (precision*recall) / (precision+recall)
perf.metrics_RF$F.measure <- 2 * precision.RF.wa * sensitivity.RF.wa / (precision.RF.wa + sensitivity.RF.wa)
perf.metrics_RF

###############################################################################
## 
## SVM for multiclass classification
## Polynomial, Radial
## 
###############################################################################

### SVC, SVM poly, SVM radial
tune_params <- c(0.1,0.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
n1 <- length(tune_params)
svc_err <- rep(NA, n1)
svm1_err <- rep(NA, n1)
svm2_err <- rep(NA, n1)

for (ii in 1:n1) {
  ### (1) support vector classifier
  svcfit <- svm(Kingdom ~ ., type="C-classification", kernel = "linear", data = data.train, cost = tune_params[ii])
  svc_pred <- predict(svcfit, newdata = data.validate)
  # Calculate the test misclassification error rate
  svc_err[ii] <- mean(svc_pred != data.validate$Kingdom)
  
  ### (2)  support vector machine with polynomial kernel of degree 3
  svmfit <- svm(Kingdom ~ ., type="C-classification", kernel = "polynomial", degree = 3, data = data.train, cost = tune_params[ii])
  svm_pred <- predict(svmfit, newdata = data.validate)
  # Calculate the test misclassification error rate
  svm1_err[ii] <- mean(svm_pred != data.validate$Kingdom)
  
  ### (3)  support vector machine with radial kernel and gamma equal to 0.1
  svmfit2 <- svm(Kingdom ~ ., type="C-classification", kernel = "radial", gamma = 0.1, data = data.train, cost = tune_params[ii])
  svm_pred2 <- predict(svmfit2, newdata = data.validate)
  # Calculate the test misclassification error rate
  svm2_err[ii] <- mean(svm_pred2 != data.validate$Kingdom)
}

plot(tune_params, svc_err, main="error v. tune params", ylab="error", ylim=c(0, 0.35))
points(svm1_err, col="red")
points(svm2_err, col="blue")

### SVMs

set.seed(228)
# Define the tuning grid for SVM with polynomial kernel
tuneGridPoly <- expand.grid(
  degree = 1:5,
  cost = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1:20)
)
# Tune the hyperparameters for SVM with polynomial kernel
svmPolyTune <- tune(
  svm,
  Kingdom ~ .,
  data = data.train,
  kernel = "polynomial",
  ranges = list(degree = 1:5, cost = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1:20)),
  tunecontrol = tune.control(sampling = "cross"),
  grid = tuneGridPoly
)

# Define the tuning grid for SVM with radial kernel
tuneGridRadial <- expand.grid(
  gamma = c(0.01,0.05, 0.1, 0.5, 1),
  cost = c(0.1, 1, 10)
)

# Tune the hyperparameters for SVM with radial kernel
svmRadialTune <- tune(
  svm,
  Kingdom ~ .,
  data = data.train,
  kernel = "radial",
  ranges = list(gamma = c(0.01,0.05, 0.1, 0.5, 1), cost = c(0.1, 1, 10)),
  #ranges = list(gamma = c(0.01,0.05, 0.1, 0.5, 1), cost = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1:20)),
  tunecontrol = tune.control(sampling = "cross"),
  grid = tuneGridRadial
)

# Train SVM with the tuned hyperparameters
svmPoly <- svm(
  Kingdom ~ .,
  data = data.train,
  kernel = "polynomial",
  degree = svmPolyTune$best.parameters$degree,
  cost = svmPolyTune$best.parameters$cost
)

svmRadial <- svm(
  Kingdom ~ .,
  data = data.train,
  kernel = "radial",
  gamma = svmRadialTune$best.parameters$gamma,
  cost = svmRadialTune$best.parameters$cost
)

# # Make predictions on the test set
predPoly <- predict(svmPoly, data.validate)
predRadial <- predict(svmRadial, data.validate)

table(predPoly, data.validate$Kingdom)
confusion.poly <- table(predPoly, data.validate$Kingdom)

table(predRadial, data.validate$Kingdom)
confusion.radial <- table(predRadial, data.validate$Kingdom)


# # Evaluate the accuracy of the models
accPoly <- sum(predPoly == data.validate$Kingdom)/nrow(data.validate)
accRadial <- sum(predRadial == data.validate$Kingdom)/nrow(data.validate)
# 

# Polynomial Kernel
### misclassification
misclassification.poly <- (confusion.poly[1,2]+ confusion.poly[1,3]+ confusion.poly[1,4]
                         + confusion.poly[2,1]+ confusion.poly[2,3]+ confusion.poly[2,4]
                         + confusion.poly[3,1]+ confusion.poly[3,2]+ confusion.poly[2,4]
                         + confusion.poly[4,1]+ confusion.poly[4,2]+ confusion.poly[4,3])/nrow(data.validate)
misclassification.poly

sensitivity.poly <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (sum(confusion.poly[,i]) == 0){
    sensitivity.poly[i] <- 0.0
  }
  else{
    sensitivity.poly[i] <- confusion.poly[,i][i]/sum(confusion.poly[,i])
  }
}

sensitivity.poly.wa <-
  (
    sensitivity.poly[1] * sum(data.validate$Kingdom == "mam") + sensitivity.poly[2] *
      sum(data.validate$Kingdom == "pri") + sensitivity.poly[3] *
      sum(data.validate$Kingdom == "rod") + sensitivity.poly[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]

# 1st element in vector is for mam, 2nd is for pri, 3rd is for rod, 4th is for vrt

### specificity = true neg / (true neg + false pos)
TN.mam <- confusion.poly[2, 2] + confusion.poly[2, 3] + confusion.poly[2, 4] + confusion.poly[3, 2] +
  confusion.poly[3, 3] + confusion.poly[3, 4] + confusion.poly[4, 2] + confusion.poly[4, 3] +
  confusion.poly[4, 4]
FP.mam <- confusion.poly[1,2]+confusion.poly[1,3]+confusion.poly[1,4]
specificity.mam <- TN.mam / (TN.mam+FP.mam)

TN.pri <-
  confusion.poly[1, 1] + confusion.poly[1, 3] + confusion.poly[1, 4] + confusion.poly[3, 1] +
  confusion.poly[3, 3] + confusion.poly[3, 4] + confusion.poly[4, 1] + confusion.poly[4, 3] +
  confusion.poly[4, 4]
FP.pri <- confusion.poly[2, 1] + confusion.poly[2, 3] + confusion.poly[2, 4] 
specificity.pri <- TN.pri / (TN.pri+FP.pri)

TN.rod <- confusion.poly[1, 1] + confusion.poly[1, 2] + confusion.poly[1, 4] + confusion.poly[2, 1] +
  confusion.poly[2, 2] + confusion.poly[2, 4] + confusion.poly[4, 1] + confusion.poly[4, 2] +
  confusion.poly[4, 4]
FP.rod <- confusion.poly[3, 1] + confusion.poly[3, 2] + confusion.poly[3, 4]
specificity.rod <- TN.rod / (TN.rod+FP.rod)

TN.vrt <-
  confusion.poly[1, 1] + confusion.poly[1, 2] + confusion.poly[1, 3] + confusion.poly[2, 1] +
  confusion.poly[2, 2] + confusion.poly[2, 3] + confusion.poly[3, 1] + confusion.poly[3, 2] + confusion.poly[3, 3]
FP.vrt <- confusion.poly[4, 1] + confusion.poly[4, 2] + confusion.poly[4, 3]
specificity.vrt <- TN.vrt / (TN.vrt+FP.vrt)
specificity.poly <- c(specificity.mam, specificity.pri, specificity.rod, specificity.vrt)

specificity.poly.wa <-
  (
    specificity.poly[1] * sum(data.validate$Kingdom == "mam") + specificity.poly[2] *
      sum(data.validate$Kingdom == "pri") + specificity.poly[3] *
      sum(data.validate$Kingdom == "rod") + specificity.poly[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]


precision.poly <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (sum(confusion.poly[i,]) == 0){
    precision.poly[i] <- 0.0
  }
  else{
    precision.poly[i] <- confusion.poly[i,][i]/sum(confusion.poly[i,])
  }
}

precision.poly.wa <-
  (
    precision.poly[1] * sum(data.validate$Kingdom == "mam") + precision.poly[2] *
      sum(data.validate$Kingdom == "pri") + precision.poly[3] *
      sum(data.validate$Kingdom == "rod") + precision.poly[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]


## F Measure
fmeasure.poly <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (precision.poly[i]+sensitivity.poly[i] == 0){
    fmeasure.poly[i] <- 0.0
  }
  else{
    fmeasure.poly[i] <- 2*(precision.poly[i]*sensitivity.poly[i])/(precision.poly[i]+sensitivity.poly[i])
  }
}
fmeasure.poly


perf.metrics_poly <- data.frame(cbind(misclassification.poly, sensitivity.poly.wa, 
                                    specificity.poly.wa, precision.poly.wa))

### F-measure = (precision*recall) / (precision+recall)
perf.metrics_poly$F.measure <- 2 * precision.poly.wa * sensitivity.poly.wa / (precision.poly.wa + sensitivity.poly.wa)
perf.metrics_poly


# Radial Kernel
### misclassification
misclassification.radial <- (confusion.radial[1,2]+ confusion.radial[1,3]+ confusion.radial[1,4]
                         + confusion.radial[2,1]+ confusion.radial[2,3]+ confusion.radial[2,4]
                         + confusion.radial[3,1]+ confusion.radial[3,2]+ confusion.radial[2,4]
                         + confusion.radial[4,1]+ confusion.radial[4,2]+ confusion.radial[4,3])/nrow(data.validate)
misclassification.radial

sensitivity.radial <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (sum(confusion.radial[,i]) == 0){
    sensitivity.radial[i] <- 0.0
  }
  else{
    sensitivity.radial[i] <- confusion.radial[,i][i]/sum(confusion.radial[,i])
  }
}

sensitivity.radial.wa <-
  (
    sensitivity.radial[1] * sum(data.validate$Kingdom == "mam") + sensitivity.radial[2] *
      sum(data.validate$Kingdom == "pri") + sensitivity.radial[3] *
      sum(data.validate$Kingdom == "rod") + sensitivity.radial[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]

# 1st element in vector is for mam, 2nd is for pri, 3rd is for rod, 4th is for vrt

### specificity = true neg / (true neg + false pos)
TN.mam <- confusion.radial[2, 2] + confusion.radial[2, 3] + confusion.radial[2, 4] + confusion.radial[3, 2] +
  confusion.radial[3, 3] + confusion.radial[3, 4] + confusion.radial[4, 2] + confusion.radial[4, 3] +
  confusion.radial[4, 4]
FP.mam <- confusion.radial[1,2]+confusion.radial[1,3]+confusion.radial[1,4]
specificity.mam <- TN.mam / (TN.mam+FP.mam)

TN.pri <-
  confusion.radial[1, 1] + confusion.radial[1, 3] + confusion.radial[1, 4] + confusion.radial[3, 1] +
  confusion.radial[3, 3] + confusion.radial[3, 4] + confusion.radial[4, 1] + confusion.radial[4, 3] +
  confusion.radial[4, 4]
FP.pri <- confusion.radial[2, 1] + confusion.radial[2, 3] + confusion.radial[2, 4] 
specificity.pri <- TN.pri / (TN.pri+FP.pri)

TN.rod <- confusion.radial[1, 1] + confusion.radial[1, 2] + confusion.radial[1, 4] + confusion.radial[2, 1] +
  confusion.radial[2, 2] + confusion.radial[2, 4] + confusion.radial[4, 1] + confusion.radial[4, 2] +
  confusion.radial[4, 4]
FP.rod <- confusion.radial[3, 1] + confusion.radial[3, 2] + confusion.radial[3, 4]
specificity.rod <- TN.rod / (TN.rod+FP.rod)

TN.vrt <-
  confusion.radial[1, 1] + confusion.radial[1, 2] + confusion.radial[1, 3] + confusion.radial[2, 1] +
  confusion.radial[2, 2] + confusion.radial[2, 3] + confusion.radial[3, 1] + confusion.radial[3, 2] + confusion.radial[3, 3]
FP.vrt <- confusion.radial[4, 1] + confusion.radial[4, 2] + confusion.radial[4, 3]
specificity.vrt <- TN.vrt / (TN.vrt+FP.vrt)
specificity.radial <- c(specificity.mam, specificity.pri, specificity.rod, specificity.vrt)

specificity.radial.wa <-
  (
    specificity.radial[1] * sum(data.validate$Kingdom == "mam") + specificity.radial[2] *
      sum(data.validate$Kingdom == "pri") + specificity.radial[3] *
      sum(data.validate$Kingdom == "rod") + specificity.radial[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]


precision.radial <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (sum(confusion.radial[i,]) == 0){
    precision.radial[i] <- 0.0
  }
  else{
    precision.radial[i] <- confusion.radial[i,][i]/sum(confusion.radial[i,])
  }
}

precision.radial.wa <-
  (
    precision.radial[1] * sum(data.validate$Kingdom == "mam") + precision.radial[2] *
      sum(data.validate$Kingdom == "pri") + precision.radial[3] *
      sum(data.validate$Kingdom == "rod") + precision.radial[4] *
      sum(data.validate$Kingdom == "vrt")
  ) / dim(data.validate)[1]


## F Measure
fmeasure.radial <- rep(0, 4)
for (i in 1:4){
  print(i)
  if (precision.radial[i]+sensitivity.radial[i] == 0){
    fmeasure.radial[i] <- 0.0
  }
  else{
    fmeasure.radial[i] <- 2*(precision.radial[i]*sensitivity.radial[i])/(precision.radial[i]+sensitivity.radial[i])
  }
}
fmeasure.radial


perf.metrics_radial <- data.frame(cbind(misclassification.radial, sensitivity.radial.wa, 
                                    specificity.radial.wa, precision.radial.wa))

### F-measure = (precision*recall) / (precision+recall)
perf.metrics_radial$F.measure <- 2 * precision.radial.wa * sensitivity.radial.wa / (precision.radial.wa + sensitivity.radial.wa)
perf.metrics_radial