library(rcompanion)
data = read.csv(file="d:\\Users\\Dell\\Desktop\\ST6090\\sarcoma.csv", header=TRUE)
dat=na.omit(data)
cramerV(dat$sex,dat$tumor.subtype)
dat$grade = NULL
dat$sex = NULL
dat$tumor.subtype = NULL
cor(dat$suv.mean, dat$suv.peak.local.suv)

# finding correlation
set.seed(4060)
options("max.print" = 1000000)
corr_of_df = cor(dat)
for (i in 1:nrow(corr_of_df)){
  correlations <- which((corr_of_df[i,] > .80) & (corr_of_df[i,] != 1))
  
  if(length(correlations) > 0){
    print(colnames(dat)[i])
    print(correlations)
  }
}
for (i in 1:nrow(corr_of_df)){
  correlations <- which((corr_of_df[i,] < -.80))
  
  if(length(correlations) > 0){
    print(colnames(dat)[i])
    print(correlations)
  }
}
#############################################
# removing correlated variables
cor_matrix <- cor(dat)


cor_matrix
# Modify correlation matrix
cor_matrix_rm <- cor_matrix
cor_matrix_rm[upper.tri(cor_matrix_rm)] <- 0
diag(cor_matrix_rm) <- 0
cor_matrix_rm
# Remove highly correlated variables
data_new <- dat[ , !apply(cor_matrix_rm,
                          2,
                          function(x) any(x > 0.95))]
head(data_new)
names(data_new)