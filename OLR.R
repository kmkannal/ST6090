library(foreign)
library(Hmisc)
library(reshape2)
library(MASS)
library(caret)
library(ggpubr)
data = read.csv(file="d:\\Users\\Dell\\Desktop\\ST6090\\sarcoma_uncor.csv", header=TRUE)
dat=na.omit(data)
dat$grade = as.factor(ifelse(dat$grade=="high",1,0))
levels(dat$grade) = c("low","high")
dat$sex=as.numeric(as.factor(dat$sex))
dat$tumor.subtype=as.numeric(as.factor(dat$tumor.subtype))
head(dat$sex)
library(dplyr)
group_by(dat, grade) %>%
  summarise(
    count = n(),
    mean = mean(grad1, na.rm = TRUE),
    sd = sd(grad1, na.rm = TRUE)
  )
dat=scale(dat)
model = glm(dat$grade~.,data = dat)
summary(model)
#describe(dat)
ggboxplot(dat, x = "grade", y = "grad1", 
          color = "grade", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
          ylab = "grad1", xlab = "grade")
# Compute the analysis of variance
res.aov <- aov(age ~ grade, data = dat)
# Summary of the analysis
summary(res.aov)

cor(dat$mean.suv,dat$mean.seg)

?p.adjust()

library(randomForest)
library(pROC)
set.seed(6041)
N = nrow(dat)
itrain = sample(1:N,N,replace=TRUE)
dat.train = dat[itrain,] 
dat.test = dat[-itrain,] 
# grow a forest:
rf.out = randomForest(dat$grade~tumor.subtype+H0+grad0.4+suv.mean+sphericity+pca.flatness+median_HIST+inv.var_GLCM+info.corr.1_GLCM+info.corr.2_GLCM, dat.train)
# fitted values for "training set"
rf.yhat = predict(rf.out, dat.train, type="class")
rf.pred = predict(rf.out, dat.test, type="class")
(tb.rf = table(rf.pred, dat.test$grade))
accuracy = sum(diag(tb.rf))/sum(tb.rf)
dat.test$grade=as.numeric(dat.test$grade)
class(dat.test$grade)
auc.rf = multiclass.roc(dat.test$grade, as.numeric(rf.pred))$auc

x.train = dat.train[,-1]
y.train = dat.train[,1]
set.seed(4061)
subsets <- c(1:10,20,30,40,ncol(dat))
ctrl <- rfeControl(functions = rfFuncs,
                   #method = "cv",
                   method = "repeatedcv",
                   repeats = 5,
                   number = 5,
                   verbose = FALSE)
rf.rfe <- rfe(x.train, y.train,
              sizes = subsets,
              rfeControl = ctrl)
predictors(rf.rfe)
varImp(rf.rfe)

library(leaps)
set.seed(4061)
data_uncor=read.csv(file="d:\\Users\\Dell\\Desktop\\ST6090\\sarcoma_uncor.csv", header=TRUE)
dat_uncor=na.omit(data_uncor)
dat_uncor$grade=as.factor(dat_uncor$grade)
dat_uncor$sex=as.numeric(as.factor(dat_uncor$sex))
dat_uncor$tumor.subtype=as.numeric(as.factor(dat_uncor$tumor.subtype))
reg.full = regsubsets(dat_uncor$grade~., data=dat_uncor, method="exhaustive",really.big=T,nvmax=10)
names(summary(reg.full))
summary(reg.full)
summary(reg.full)$which #tracks covariate selection in 8 subsets
summary(reg.full)$outmat
?regsubsets
set.seed(4061)
K=10
final=matrix(nrow=K,106)
for(k in 1:K){
  reg.fwd = regsubsets(dat$grade~., data=dat,method="forward")
  isel = which(summary(reg.fwd)$outmat=="*")
  final[k,isel] = 1
  final[k,-isel] = 0
}
names(summary(reg.fwd))$outmat
set.seed(4061)
reg.bwd = regsubsets(dat_uncor$grade~., data=dat_uncor, nvmax=10,method="backward")
summary(reg.bwd)$outmat
