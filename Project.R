library(caret)
library(randomForest)
library(e1071)
library(class)
Sarcoma = read.csv(file="d:\\Users\\Dell\\Desktop\\ST6090\\sarcoma.csv", header=TRUE)
Sarcoma$grade = as.factor(ifelse(Sarcoma$grade=="high",1,0))
levels(Sarcoma$grade) = c("low","high")
Sar=na.omit(Sarcoma)
Sar$sex=as.numeric(as.factor(Sar$sex))
Sar$tumor.subtype=as.numeric(as.factor(Sar$tumor.subtype))

summary(Sar$surface)

#Repeated 5 fold CV
train_control <- trainControl(method = "repeatedcv",
                              number = 5, repeats = 3)
model <- train(grade~., data = Sar,
               trControl = train_control, method = "rf")

print(model)

model_1=train(grade~H0+H1+grad0.25+grad0.3+grad0.4+grad0.5+grad0.6+grad1+suv.mean+tlg.seg+suv.peak.local.seg+
  volume+surface+sv.ratio+sphericity+asphericity+major.axis.length+minor.axis.length+
  least.axis.length+pca.elongation+pca.flatness+energy+skew_HIST+kurt_HIST+median_HIST+p10_HIST+
  p90_HIST+mode_HIST+medad_HIST+CoV_HIST+qcod_HIST+uniformity_HIST+max.grad_HIST+
  Reg.max.grad.grey_HIST+min.grad_HIST+Reg.min.grad.grey_HIST+entropy_GLCM+diff.var_GLCM+
  sum.var_GLCM+sum.entropy_GLCM+energy_GLCM+dissimilarity_GLCM+inv.diff.mom_GLCM+
  inv.diff.mom.norm_GLCM+inv.var_GLCM+correlation_GLCM+auto.corr_GLCM+clust.shade_GLCM+
  clust.prom_GLCM+info.corr.1_GLCM+info.corr.2_GLCM+Reg.Het1+Reg.Het0+Reg.grad.min+
  Reg.grad.0.05+Reg.grad.0.95+Reg.grad.1+age+sex+tumor.subtype, data = Sar,
  trControl = train_control, method = "rf")
print(model_1)
model_2=train(grade~suv.mean+suv.peak.local.suv+max.suv+
                mean.seg+suv.peak.local.seg+tumor.subtype+       
                tlg.seg+age+grad1+max.grad_HIST+grad0.6+Reg.grad.min+
                inv.var_GLCM+pca.flatness+H1+
                volume+grad0.95+pca.elongation+
                grad0.8+info.corr.2_GLCM+grad0.5+
                diff.var_GLCM+grad0.7+grad0.85+
                grad0.9+grad0.75+joint.var_GLCM+
                dissimilarity_GLCM+Reg.grad.0.05+homogeneity.norm_GLCM,data = Sar,
              trControl = train_control, method = "rf")
print(model_2)
model_3=train(grade~suv.mean+suv.peak.local.seg+tumor.subtype+tlg.seg+age+grad0.6+grad1+inv.var_GLCM,data = Sar,
              trControl = train_control, method = "rf")
print(model_3)

set.seed(4061)
n = nrow(Sar)
dat = Sar[sample(1:n, n, replace=FALSE), ] #shuffle the dataset

# get a random training sample containing 70% of original sample:
i.cv = sample(1:n, round(.7*n), replace=FALSE)
dat.cv = dat[i.cv,] # use this for CV (train+test)
dat.valid = dat[-i.cv,] #hold out sample
# perform repeated K-fold CV:
R=3
r=0
K = 5 #folds
K.knn = 3
N = length(i.cv)
folds = cut(1:N, K, labels=FALSE)
#accuracy
acc.glm = acc.svm = acc.knn = acc.rf = numeric(K)
pres.glm = pres.svm = pres.knn = pres.rf = numeric(K)
sens.glm = sens.svm = sens.knn = sens.rf = numeric(K)
f1.glm = f1.svm = f1.knn = f1.rf = numeric(K)
while (r <= R){
  for(k in 1:K){     
    #split into train and test samples:
    i.train	= which(folds!=k)
    dat.train = dat.cv[i.train,]
    dat.test = dat.cv[-i.train,]
    x.train = dat.train[,c("tumor.subtype","H1","max.suv","mean.seg","suv.mean","suv.peak.local.seg","age","tlg.seg","grad0.5","grad0.6","Reg.max.grad.grey_HIST","Reg.grad.1","grad0.7","pca.elongation","grad0.75","volume","correlation_GLCM","pca.flatness","Reg.grad.min","max.grad_HIST","grad0.95","grad1","grad0.85","grad0.8","info.corr.2_GLCM","Reg.grad.0.05","qcod_HIST","CoV_HIST")]
    y.train = dat.train[,1]
    x.test = dat.test[,c("tumor.subtype","H1","max.suv","mean.seg","suv.mean","suv.peak.local.seg","age","tlg.seg","grad0.5","grad0.6","Reg.max.grad.grey_HIST","Reg.grad.1","grad0.7","pca.elongation","grad0.75","volume","correlation_GLCM","pca.flatness","Reg.grad.min","max.grad_HIST","grad0.95","grad1","grad0.85","grad0.8","info.corr.2_GLCM","Reg.grad.0.05","qcod_HIST","CoV_HIST")]
    y.test = dat.test[,1]
    #GLM
    glm.o = glm(grade~tumor.subtype+H1+max.suv+mean.seg+suv.mean+suv.peak.local.seg+age+tlg.seg+grad0.5+grad0.6+Reg.max.grad.grey_HIST+Reg.grad.1+grad0.7+pca.elongation+grad0.75+volume+correlation_GLCM+pca.flatness+Reg.grad.min+max.grad_HIST+grad0.95+grad1+grad0.85+grad0.8+info.corr.2_GLCM+Reg.grad.0.05+qcod_HIST+CoV_HIST, data=dat.train, family=binomial(logit))
    glm.p = (predict(glm.o, newdata=dat.test, type="response")>0.5)
    tb.glm = table(glm.p, y.test)
    acc.glm[k] = sum(diag(tb.glm)) / sum(tb.glm)
    pres.glm[k] = tb.glm[2,2]/sum(tb.glm[2,])
    sens.glm[k] = tb.glm[2,2]/sum(tb.glm[,2])
    f1.glm[k] = (2*pres.glm[k]*sens.glm[k])/(pres.glm[k]+sens.glm[k])
    #KNN
    knn.o = knn(x.train, x.test, y.train, K.knn)
    knn.p = knn.o
    tb.knn = table(knn.p, y.test)
    acc.knn[k] = sum(diag(tb.knn)) / sum(tb.knn)
    pres.knn[k] = tb.knn[2,2]/sum(tb.knn[2,])
    sens.knn[k] = tb.knn[2,2]/sum(tb.knn[,2]) #recall
    f1.knn[k] = (2*pres.knn[k]*sens.knn[k])/(pres.knn[k]+sens.knn[k])
    #SVM
    svmo.lin = svm(grade~tumor.subtype+H1+max.suv+mean.seg+suv.mean+suv.peak.local.seg+age+tlg.seg+grad0.5+grad0.6+Reg.max.grad.grey_HIST+Reg.grad.1+grad0.7+pca.elongation+grad0.75+volume+correlation_GLCM+pca.flatness+Reg.grad.min+max.grad_HIST+grad0.95+grad1+grad0.85+grad0.8+info.corr.2_GLCM+Reg.grad.0.05+qcod_HIST+CoV_HIST, data=dat.train, kernel='radial')
    pred.lin = predict(svmo.lin, newdata=dat.test)
    #Confusion matrix
    #tb.svm = table(svmy.lin, y.test)
    tb.svm = table(pred.lin, y.test)
    acc.svm[k] = sum(diag(tb.svm)) / sum(tb.svm)
    pres.svm[k] = tb.svm[2,2]/sum(tb.svm[2,])
    sens.svm[k] = tb.svm[2,2]/sum(tb.svm[,2]) 
    f1.svm[k] = (2*pres.svm[k]*sens.svm[k])/(pres.svm[k]+sens.svm[k])
    #RF
    rf.out = randomForest(grade~tumor.subtype+H1+max.suv+mean.seg+suv.mean+suv.peak.local.seg+age+tlg.seg+grad0.5+grad0.6+Reg.max.grad.grey_HIST+Reg.grad.1+grad0.7+pca.elongation+grad0.75+volume+correlation_GLCM+pca.flatness+Reg.grad.min+max.grad_HIST+grad0.95+grad1+grad0.85+grad0.8+info.corr.2_GLCM+Reg.grad.0.05+qcod_HIST+CoV_HIST, dat.train)
    rf.yhat = predict(rf.out, dat.train, type="class")
    # fitted values for "test set"
    rf.pred = predict(rf.out, dat.test, type="class")
    tb.rf = table(rf.pred, y.test)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    pres.rf[k] = tb.rf[2,2]/sum(tb.rf[2,])
    sens.rf[k] = tb.rf[2,2]/sum(tb.rf[,2])
    f1.rf[k] = (2*pres.rf[k]*sens.rf[k])/(pres.rf[k]+sens.rf[k])
    }
  r <- r+1
}
boxplot(as.vector(acc.svm), as.vector(acc.glm), as.vector(acc.knn),
        as.vector(acc.rf), names=c("SVM","GLM","KNN","RF"), 
        main="Accuracy for Backward Stepwise")
round(mean(acc.glm),2)
round(mean(acc.knn),2)
round(mean(acc.svm),2)
round(mean(acc.rf),2)
round(sd(acc.glm),2)
round(sd(acc.knn),2)
round(sd(acc.svm),2)
round(sd(acc.rf),2)

round(0.6818182,2)

set.seed(4061)
rf.out = randomForest(grade~., dat.train)
varImpPlot(rf.out,pch=15,main="Important variables")

library(OneR)
library(FSelectorRcpp)
library(praznik)

model <- OneR(Sar, verbose = TRUE)
oner_Type <- OneR::OneR(Sar$grade~Sar$tumor.subtype)
oner_Type

infogain_Type <- information_gain(grade~., data=dat.train, type="infogain")
infogain_Type

require(dplyr)
df1 %>% count(Sarcoma$grade)

boxplot(Sar$grade,Sar$pca.flatness, xlab = "Grade",
        ylab = "pca.flatness",names=c("Low","High"))
levels(Sar$grade)

set.seed(4061)
n = nrow(Sar)
dat = Sar[sample(1:n, n, replace=FALSE), ] #shuffle the dataset

# get a random training sample containing 70% of original sample:
i.cv = sample(1:n, round(.7*n), replace=FALSE)
dat.cv = dat[i.cv,] # use this for CV (train+test)
dat.valid = dat[-i.cv,] #hold out sample
# perform repeated K-fold CV:
R=3
r=0
K = 5 #folds
K.knn = 3
N = length(i.cv)
folds = cut(1:N, K, labels=FALSE)
#accuracy
acc.rf = acc.rf.pre = acc.rf.rfe = acc.fin = numeric(K)
while (r <= R){
  for(k in 1:K){     
    #split into train and test samples:
    i.train	= which(folds!=k)
    dat.train = dat.cv[i.train,]
    dat.test = dat.cv[-i.train,]
    y.test = dat.test[,1]
    #RF
    rf.out = randomForest(grade~., dat.train)
    rf.yhat = predict(rf.out, dat.train, type="class")
    # fitted values for "test set"
    rf.pred = predict(rf.out, dat.test, type="class")
    tb.rf = table(rf.pred, y.test)
    acc.rf[k] = sum(diag(tb.rf)) / sum(tb.rf)
    #Pre-filter
    rf.out.pre = randomForest(grade~H0+H1+grad0.25+grad0.3+grad0.4+grad0.5+grad0.6+grad1+suv.mean+tlg.seg+suv.peak.local.seg+volume+surface+sv.ratio+sphericity+asphericity+major.axis.length+minor.axis.length+least.axis.length+pca.elongation+pca.flatness+energy+skew_HIST+kurt_HIST+median_HIST+p10_HIST+p90_HIST+mode_HIST+medad_HIST+CoV_HIST+qcod_HIST+uniformity_HIST+max.grad_HIST+Reg.max.grad.grey_HIST+min.grad_HIST+Reg.min.grad.grey_HIST+entropy_GLCM+diff.var_GLCM+sum.var_GLCM+sum.entropy_GLCM+energy_GLCM+dissimilarity_GLCM+inv.diff.mom_GLCM+inv.diff.mom.norm_GLCM+inv.var_GLCM+correlation_GLCM+auto.corr_GLCM+clust.shade_GLCM+clust.prom_GLCM+info.corr.1_GLCM+info.corr.2_GLCM+Reg.Het1+Reg.Het0+Reg.grad.min+Reg.grad.0.05+Reg.grad.0.95+Reg.grad.1+age+sex+tumor.subtype, dat.train)
    rf.yhat.pre = predict(rf.out.pre, dat.train, type="class")
    # fitted values for "test set"
    rf.pred.pre = predict(rf.out.pre, dat.test, type="class")
    tb.rf.pre = table(rf.pred.pre, y.test)
    acc.rf.pre[k] = sum(diag(tb.rf.pre)) / sum(tb.rf.pre)
  }
  r <- r+1
}
mean(acc.rf)
median(acc.rf)
boxplot(as.vector(acc.rf),as.vector(acc.rf.pre), names=c("Whole","Prefiltering"), 
        main="Accuracy Boxplot")

