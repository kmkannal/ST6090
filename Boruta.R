library(Boruta)
library(caret)
sarcoma.grade=read.csv(file="d:\\Users\\Dell\\Desktop\\ST6090\\sarcoma_grade.csv", header=TRUE)
sarcoma.grade$sex = as.numeric(as.factor(sarcoma.grade$sex))
sarcoma.grade$grade = as.numeric(as.factor(sarcoma.grade$grade))
sarcoma.grade.cleaned=na.omit(sarcoma.grade)
scale(sarcoma.grade.cleaned, center = TRUE, scale = TRUE)
set.seed(4060)
p=ncol(sarcoma.grade.cleaned)
n=nrow(sarcoma.grade.cleaned)
i.train = sample(1:n, size=108, replace=FALSE)
dat.train = sarcoma.grade.cleaned[i.train,]
dat.test = sarcoma.grade.cleaned[-i.train,]
K = 10
folds = cut(1:n, breaks=K, labels=FALSE)
final.boruta=matrix(nrow=K,ncol=106)
final.boruta
for(k in 1:K){
  itr = which(folds!=k)
  boruta.train <- Boruta(grade~., data=dat.train, doTrace = 2)
  isel = which(boruta.train$finalDecision=="Confirmed")
  final.boruta[k,isel] = 1
  final.boruta[k,-isel] = 0
  #features[k]=getSelectedAttributes(final.boruta, withTentative = F)
}
#rpart
library(caret)
set.seed(100)
rPartgrade <- train(grade~., data=sarcoma.grade.cleaned, method="rpart")
rpartgradeImp <- varImp(rPartgrade)
print(rpartgradeImp)
#lasso
library(glmnet)
x <- as.matrix(sarcoma.grade.cleaned[,-1])
y <- as.matrix(sarcoma.grade.cleaned[,1])
set.seed(100)
cv.grade.lasso <- cv.glmnet(x, y, family='multinomial',nfolds = 5)
cv.grade=glmnet(x, y, family='multinomial', alpha=1, lambda = cv.grade.lasso$lambda.min)
coef(cv.grade)
#tumor subtype as response variable
sarcoma.subtp=read.csv(file="d:\\Users\\Dell\\Desktop\\ST6090\\sarcoma_subtp.csv", header=TRUE)
sarcoma.subtp$sex = as.numeric(as.factor(sarcoma.subtp$sex))
sarcoma.subtp$tumor.subtype = as.numeric(as.factor(sarcoma.subtp$tumor.subtype))
sarcoma.subtp.cleaned=na.omit(sarcoma.subtp)
boruta.train1 <- Boruta(tumor.subtype~., data=sarcoma.subtp.cleaned, doTrace = 2)
print(boruta.train1)
final.boruta1 <- TentativeRoughFix(boruta.train1)
print(final.boruta1)
getSelectedAttributes(final.boruta1, withTentative = F)
#lasso
x <- as.matrix(sarcoma.subtp.cleaned[,-1])
y <- as.matrix(sarcoma.subtp.cleaned[,1])
set.seed(100)
cv.subtp.lasso <- cv.glmnet(x, y, family='multinomial',nfolds = 5)
cv.subtp=glmnet(x, y, family='multinomial', alpha=1, lambda = cv.subtp.lasso$lambda.min)
coef(cv.subtp)
#surtim as response variable
sarcoma.surtim=read.csv(file="d:\\Users\\Dell\\Desktop\\ST6090\\sarcoma_surtim.csv", header=TRUE)
sarcoma.surtim$sex = as.numeric(as.factor(sarcoma.surtim$sex))
sarcoma.surtim.cleaned=na.omit(sarcoma.surtim)
boruta.train2 <- Boruta(surtim~., data=sarcoma.surtim.cleaned, doTrace = 2)
print(boruta.train2)
final.boruta2 <- TentativeRoughFix(boruta.train2)
print(final.boruta2)
getSelectedAttributes(final.boruta2, withTentative = F)
#surind as response variable
sarcoma.surind=read.csv(file="d:\\Users\\Dell\\Desktop\\ST6090\\sarcoma_surind.csv", header=TRUE)
sarcoma.surind$sex = as.numeric(as.factor(sarcoma.surind$sex))
sarcoma.surind.cleaned=na.omit(sarcoma.surind)
boruta.train3 <- Boruta(surind~., data=sarcoma.surind.cleaned, doTrace = 2)
print(boruta.train3)
final.boruta3 <- TentativeRoughFix(boruta.train3)
print(final.boruta3)
getSelectedAttributes(final.boruta3, withTentative = F)
#lasso
x <- as.matrix(sarcoma.surind.cleaned[,-1])
y <- as.matrix(sarcoma.surind.cleaned[,1])
set.seed(100)
cv.surind.lasso <- cv.glmnet(x, y, family='binomial',nfolds = 5)
cv.surind=glmnet(x, y, family='binomial', alpha=1, lambda = cv.surind.lasso$lambda.min)
coef(cv.surind)

#Univariate analysis
fit4 <- manova(cbind(sex
                     ,H0
                     ,H1
                     ,grad0.25
                     ,grad0.3
                     ,grad0.4
                     ,grad0.5
                     ,grad0.6
                     ,grad0.7
                     ,grad0.75
                     ,grad0.8
                     ,grad0.85
                     ,grad0.9
                     ,grad0.95
                     ,grad1
                     ,suv.mean
                     ,max.suv
                     ,suv.peak.local.suv
                     ,mean.seg
                     ,tlg.seg
                     ,suv.peak.local.seg
                     ,volume
                     ,surface
                     ,sv.ratio
                     ,compactness1
                     ,compactness2
                     ,sphericity
                     ,asphericity
                     ,l.major
                     ,l.minor
                     ,l.least
                     ,major.axis.length
                     ,minor.axis.length
                     ,least.axis.length
                     ,pca.elongation
                     ,pca.flatness
                     ,energy
                     ,RMS
                     ,mean_HIST
                     ,var_HIST
                     ,skew_HIST
                     ,kurt_HIST
                     ,median_HIST
                     ,p10_HIST
                     ,p90_HIST
                     ,mode_HIST
                     ,iqr_HIST
                     ,mad_HIST
                     ,rmad_HIST
                     ,medad_HIST
                     ,CoV_HIST
                     ,qcod_HIST
                     ,entropy_HIST
                     ,uniformity_HIST
                     ,max.grad_HIST
                     ,Reg.max.grad.grey_HIST
                     ,min.grad_HIST
                     ,Reg.min.grad.grey_HIST
                     ,joint.max_GLCM
                     ,joint.var_GLCM
                     ,entropy_GLCM
                     ,diff.var_GLCM
                     ,diff.entropy_GLCM
                     ,sum.avg_GLCM
                     ,sum.var_GLCM
                     ,sum.entropy_GLCM
                     ,energy_GLCM
                     ,contrast_GLCM
                     ,dissimilarity_GLCM
                     ,homogeneity_GLCM
                     ,homogeneity.norm_GLCM
                     ,inv.diff.mom_GLCM
                     ,inv.diff.mom.norm_GLCM
                     ,inv.var_GLCM
                     ,correlation_GLCM
                     ,auto.corr_GLCM
                     ,clust.shade_GLCM
                     ,clust.prom_GLCM
                     ,info.corr.1_GLCM
                     ,info.corr.2_GLCM
                     ,Reg.Het1
                     ,Reg.Het0
                     ,Reg.grad.min
                     ,Reg.grad.0.05
                     ,Reg.grad.0.1
                     ,Reg.grad.0.15
                     ,Reg.grad.0.2
                     ,Reg.grad.0.25
                     ,Reg.grad.0.3
                     ,Reg.grad.0.35
                     ,Reg.grad.0.4
                     ,Reg.grad.0.45
                     ,Reg.grad.0.5
                     ,Reg.grad.0.55
                     ,Reg.grad.0.6
                     ,Reg.grad.0.65
                     ,Reg.grad.0.7
                     ,Reg.grad.0.75
                     ,Reg.grad.0.8
                     ,Reg.grad.0.85
                     ,Reg.grad.0.9
                     ,Reg.grad.0.95
                     ,Reg.grad.1
                     ,age)~surind, data=sarcoma.surind.cleaned)
summary(fit4,tol=0)

install.packages("olsrr")
library(olsrr)
model1 <- lm(grade~., data = sarcoma.grade.cleaned)
ols_vif_tol(model1)
model2 <- lm(tumor.subtype~., data = sarcoma.subtp.cleaned)
ols_vif_tol(model2)
model3 <- lm(surtim~., data = sarcoma.surtim.cleaned)
ols_vif_tol(model3)
model4 <- lm(surind~., data = sarcoma.surind.cleaned)
ols_vif_tol(model4)

library(Boruta)
library(caret)
sarcoma=read.csv(file="d:\\Users\\Dell\\Desktop\\ST6090\\sarcoma.csv", header=TRUE)
sarcoma$sex = as.numeric(as.factor(sarcoma$sex))
sarcoma$grade = as.numeric(as.factor(sarcoma$grade))
sarcoma$tumor.subtype = as.numeric(as.factor(sarcoma$tumor.subtype))
sarcoma.cleaned=na.omit(sarcoma)
scale(sarcoma.cleaned, center = TRUE, scale = TRUE)
set.seed(4060)
p=ncol(sarcoma.cleaned)
n=nrow(sarcoma.cleaned)
K = 10
folds = cut(1:n, breaks=K, labels=FALSE)
final.boruta=matrix(nrow=K,ncol=p)
for(k in 1:K){
  i.train	= which(folds!=k)
  dat.train = sarcoma.cleaned[i.train, ]
  dat.test = sarcoma.cleaned[-i.train, ]
  boruta.train <- Boruta(grade~., data=dat.train, doTrace = 2)
  isel = which(boruta.train$finalDecision=="Confirmed")
  final.boruta[k,isel] = 1
  final.boruta[k,-isel] = 0
  #features[k]=getSelectedAttributes(final.boruta, withTentative = F)
}

library(glmnet)
xm = model.matrix(grade~.,data=sarcoma.cleaned)[,-1]
p
set.seed(4061) # for reproducibility
K = 10
folds = cut(1:n, breaks=K, labels=FALSE)
sel = matrix(nrow=K, ncol=105)
nms = colnames(subset(sarcoma.cleaned, select = -c(grade)))
colnames(sel)=nms
for(k in 1:K){
  itr = which(folds!=k)
  lasso.cv = cv.glmnet(xm[itr,], sarcoma.cleaned$grade[itr])
  lasso = glmnet(xm[itr,], sarcoma.cleaned$grade[itr], lambda=lasso.cv$lambda.min)
  isel = which(coef(lasso)[-1] != 0)
  apply(sel,2,mean)*100
  sel[k,isel] = 1
  sel[k,-isel] = 0
}

set.seed(6041)
#RFE 
control <- rfeControl(functions = rfFuncs, # random forest
                      method = "repeatedcv", # repeated cv
                      repeats = 5, # number of repeats
                      number = 10) # number of folds

#rfe code
inTrain <- createDataPartition(sarcoma.cleaned$grade, p = .80, list = FALSE)[,1]
x=sarcoma.cleaned[,-1]
x_train <- x[inTrain,]
x_test <- x[-inTrain,]
y_train <- sarcoma.cleaned$grade[inTrain]
y_test <- sarcoma.cleaned$grade[-inTrain]
subsets <- c(1:5, 10, 15, 20, ncol(x))
result_rfe1 <- rfe(x = x_train,y = y_train,sizes = subsets,rfeControl = control,)
predictors(result_rfe1)



