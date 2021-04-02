setwd("E:/Rdata/")
traincsv <- read.csv("E:/Rdata/train.csv")
testcsv <- read.csv("E:/Rdata/test.csv")
library(e1071);library(pROC);library(ada);library(randomForest);library(nnet);library(rpart);library(party)
trainset	<- read.csv("E:/Rdata/train.csv");
train		<- trainset[c(1:7000),]
test		<- trainset[c(7001:10000),-1]
train[,2]	<- ifelse(train[,2]>=50,mean(train[,2],trim=0.7),train[,2])
test[,1]	<- ifelse(test[,1]>=50,mean(test[,1],trim=0.7),test[,1])
train[,4]	<- ifelse(train[,4]>=95,mean(train[,4],trim=0.7),train[,4])
test[,3]	<- ifelse(test[,3]>=95,mean(test[,3],trim=0.7),test[,3])
train[,5]	<- ifelse(train[,5]>=9,mean(train[,5],trim=0.7),train[,5])
test[,4]	<- ifelse(test[,4]>=9,mean(test[,4],trim=0.7),test[,4])
train[,6]	<- ifelse(train[,6]>=10000,mean(train[,6],trim=0.7),train[,6])
test[,5]	<- ifelse(test[,5]>=10000,mean(test[,5],trim=0.7),test[,5])
train[,7]	<- ifelse(train[,7]>=20,mean(train[,7],trim=0.7),train[,7])
test[,6]	<- ifelse(test[,6]>=20,mean(test[,6],trim=0.7),test[,6])
train[,8]	<- ifelse(train[,8]>=95,mean(train[,8],trim=0.7),train[,8])
test[,7]	<- ifelse(test[,7]>=95,mean(test[,7],trim=0.7),test[,7])
train[,9]	<- ifelse(train[,9]>=5,mean(train[,9],trim=0.7),train[,9])
test[,8]	<- ifelse(test[,8]>=5,mean(test[,8],trim=0.7),test[,8])
train[,10]	<- ifelse(train[,10]>=95,mean(train[,10],trim=0.7),train[,10])
test[,9]	<- ifelse(test[,9]>=95,mean(test[,9],trim=0.7),test[,9])
train$DepRisk	<- abs(train$age-100)*train$NumberOfDependents
test$DepRisk	<- abs(test$age-100)*test$NumberOfDependents
train$DepIncRisk	<- train$MonthlyIncome*train$NumberOfDependents
test$DepIncRisk	<- test$MonthlyIncome*test$NumberOfDependents
train$InvDebtRisk <- abs(train$age-100)*train$DebtRatio
test$InvDebtRisk	<- abs(test$age-100)*test$DebtRatio
train$TDL		<- train[,4]*30+train[,8]*90+train[,10]*60
test$TDL		<- test[,3]*30+test[,7]*90+test[,9]*60
train$TDLRUOUL	<- train$TDL*train$RevolvingUtilizationOfUnsecuredLines
test$TDLRUOUL	<- test$TDL*test$RevolvingUtilizationOfUnsecuredLines
train$Money		<- (1-train$RevolvingUtilizationOfUnsecuredLines)*train$MonthlyIncome
test$Money		<- (1-test$RevolvingUtilizationOfUnsecuredLines)*test$MonthlyIncome
train$DebtInc	<- train$DebtRatio*train$MonthlyIncome
test$DebtInc	<- test$DebtRatio*test$MonthlyIncome
train$ageappdebt	<- ifelse(train$age>=90,(mean(train$DebtRatio[train$age>=90])),0) +
   ifelse((train$age<90 & train$age>=80),(mean(train$DebtRatio[train$age<90 & train$age>=80])),0) +
   ifelse((train$age<80 & train$age>=70),(mean(train$DebtRatio[train$age<80 & train$age>=70])),0) +
   ifelse((train$age<70 & train$age>=60),(mean(train$DebtRatio[train$age<70 & train$age>=60])),0) +
   ifelse((train$age<60 & train$age>=50),(mean(train$DebtRatio[train$age<60 & train$age>=50])),0) +
   ifelse((train$age<50 & train$age>=40),(mean(train$DebtRatio[train$age<50 & train$age>=40])),0) +
   ifelse((train$age<40 & train$age>=30),(mean(train$DebtRatio[train$age<40 & train$age>=30])),0) +
   ifelse((train$age<30 & train$age>=20),(mean(train$DebtRatio[train$age<30 & train$age>=20])),0)
test$ageappdebt	<- ifelse(test$age>=90,(mean(test$DebtRatio[test$age>=90])),0) +
   ifelse((test$age<90 & test$age>=80),(mean(test$DebtRatio[test$age<90 & test$age>=80])),0) +
   ifelse((test$age<80 & test$age>=70),(mean(test$DebtRatio[test$age<80 & test$age>=70])),0) +
   ifelse((test$age<70 & test$age>=60),(mean(test$DebtRatio[test$age<70 & test$age>=60])),0) +
   ifelse((test$age<60 & test$age>=50),(mean(test$DebtRatio[test$age<60 & test$age>=50])),0) +
   ifelse((test$age<50 & test$age>=40),(mean(test$DebtRatio[test$age<50 & test$age>=40])),0) +
   ifelse((test$age<40 & test$age>=30),(mean(test$DebtRatio[test$age<40 & test$age>=30])),0) +
   ifelse((test$age<30 & test$age>=20),(mean(test$DebtRatio[test$age<30 & test$age>=20])),0)
train$DebtDev	<- train$DebtRatio-train$ageappdebt
test$DebtDev	<- test$DebtRatio-test$ageappdebt
#______RANDOM FOREST: CLASS_________________________ 
set.seed(10)
RFC <- randomForest(train[,-c(1,17,19)],factor(train[,1]),mtry=2,do.trace=FALSE,importance=TRUE,ntree=500)
a <- mean(as.numeric(predict(RFC,test,type="prob")[,2]>.5)!=trainset[c(7001:10000),1]);a;round((a-0.2203333),6)
#Expected Misclass Rate: Real data=0.2203333, OOB error est 22.94%

#__________FINAL MODEL_____________________
train	<- read.csv("E:/Rdata/train.csv")
test <- read.csv("E:/Rdata/test.csv")
train[,2]	<- ifelse(train[,2]>=50,mean(train[,2],trim=0.7),train[,2])
test[,1]	<- ifelse(test[,1]>=50,mean(test[,1],trim=0.7),test[,1])
train[,4]	<- ifelse(train[,4]>=95,mean(train[,4],trim=0.7),train[,4])
test[,3]	<- ifelse(test[,3]>=95,mean(test[,3],trim=0.7),test[,3])
train[,5]	<- ifelse(train[,5]>=9,mean(train[,5],trim=0.7),train[,5])
test[,4]	<- ifelse(test[,4]>=9,mean(test[,4],trim=0.7),test[,4])
train[,6]	<- ifelse(train[,6]>=10000,mean(train[,6],trim=0.7),train[,6])
test[,5]	<- ifelse(test[,5]>=10000,mean(test[,5],trim=0.7),test[,5])
train[,7]	<- ifelse(train[,7]>=20,mean(train[,7],trim=0.7),train[,7])
test[,6]	<- ifelse(test[,6]>=20,mean(test[,6],trim=0.7),test[,6])
train[,8]	<- ifelse(train[,8]>=95,mean(train[,8],trim=0.7),train[,8])
test[,7]	<- ifelse(test[,7]>=95,mean(test[,7],trim=0.7),test[,7])
train[,9]	<- ifelse(train[,9]>=5,mean(train[,9],trim=0.7),train[,9])
test[,8]	<- ifelse(test[,8]>=5,mean(test[,8],trim=0.7),test[,8])
train[,10]	<- ifelse(train[,10]>=95,mean(train[,10],trim=0.7),train[,10])
test[,9]	<- ifelse(test[,9]>=95,mean(test[,9],trim=0.7),test[,9])
train$DepRisk	<- abs(train$age-100)*train$NumberOfDependents
test$DepRisk	<- abs(test$age-100)*test$NumberOfDependents
train$DepIncRisk	<- train$MonthlyIncome*train$NumberOfDependents
test$DepIncRisk	<- test$MonthlyIncome*test$NumberOfDependents
train$InvDebtRisk <- abs(train$age-100)*train$DebtRatio
test$InvDebtRisk	<- abs(test$age-100)*test$DebtRatio
train$TDL		<- train[,4]*30+train[,8]*90+train[,10]*60
test$TDL		<- test[,3]*30+test[,7]*90+test[,9]*60
train$TDLRUOUL	<- train$TDL*train$RevolvingUtilizationOfUnsecuredLines
test$TDLRUOUL	<- test$TDL*test$RevolvingUtilizationOfUnsecuredLines
train$Money		<- (1-train$RevolvingUtilizationOfUnsecuredLines)*train$MonthlyIncome
test$Money		<- (1-test$RevolvingUtilizationOfUnsecuredLines)*test$MonthlyIncome
train$DebtInc	<- train$DebtRatio*train$MonthlyIncome
test$DebtInc	<- test$DebtRatio*test$MonthlyIncome
train$ageappdebt	<- ifelse(train$age>=90,(mean(train$DebtRatio[train$age>=90])),0) +
   ifelse((train$age<90 & train$age>=80),(mean(train$DebtRatio[train$age<90 & train$age>=80])),0) +
   ifelse((train$age<80 & train$age>=70),(mean(train$DebtRatio[train$age<80 & train$age>=70])),0) +
   ifelse((train$age<70 & train$age>=60),(mean(train$DebtRatio[train$age<70 & train$age>=60])),0) +
   ifelse((train$age<60 & train$age>=50),(mean(train$DebtRatio[train$age<60 & train$age>=50])),0) +
   ifelse((train$age<50 & train$age>=40),(mean(train$DebtRatio[train$age<50 & train$age>=40])),0) +
   ifelse((train$age<40 & train$age>=30),(mean(train$DebtRatio[train$age<40 & train$age>=30])),0) +
   ifelse((train$age<30 & train$age>=20),(mean(train$DebtRatio[train$age<30 & train$age>=20])),0)
test$ageappdebt	<- ifelse(test$age>=90,(mean(test$DebtRatio[test$age>=90])),0) +
   ifelse((test$age<90 & test$age>=80),(mean(test$DebtRatio[test$age<90 & test$age>=80])),0) +
   ifelse((test$age<80 & test$age>=70),(mean(test$DebtRatio[test$age<80 & test$age>=70])),0) +
   ifelse((test$age<70 & test$age>=60),(mean(test$DebtRatio[test$age<70 & test$age>=60])),0) +
   ifelse((test$age<60 & test$age>=50),(mean(test$DebtRatio[test$age<60 & test$age>=50])),0) +
   ifelse((test$age<50 & test$age>=40),(mean(test$DebtRatio[test$age<50 & test$age>=40])),0) +
   ifelse((test$age<40 & test$age>=30),(mean(test$DebtRatio[test$age<40 & test$age>=30])),0) +
   ifelse((test$age<30 & test$age>=20),(mean(test$DebtRatio[test$age<30 & test$age>=20])),0)
train$DebtDev	<- train$DebtRatio-train$ageappdebt
test$DebtDev	<- test$DebtRatio-test$ageappdebt
set.seed(10)
RFC <- randomForest(train[,-c(1,17,19)],factor(train[,1]),mtry=2,do.trace=FALSE,importance=TRUE,ntree=500);RFC
preds <- as.numeric(predict(RFC,test,type="prob")[,2]>.5)

#Expected Misclass Rate: Real data=(NA,no real data), OOB error est 22.84%

dump("preds","predict.R")
