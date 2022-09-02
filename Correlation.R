# Eric Wolsztynski, UCC, May 2022
# This short script provides an example of how to keep track of selection rates from a filter in a resampling framework.
# For this example, we use bootstrapping of the correlation filter.
# Note that we transform the feature matrix using model.matrix, which recodes categorical variables, and as a result further work is needed to "recombine" the selection rates for an overall categorical variable based on the selection rates of its levels.
# Other options where no recoding happens are also "allowed"! In such cases, the original variable names would likely be preserved.

library(caret)
library(ISLR)

# load some data
data = read.csv(file="d:\\Users\\Dell\\Desktop\\ST6090\\sarcoma.csv", header=TRUE)
dat=na.omit(data)
dat$grade = NULL
dat$sex = as.numeric(as.factor(dat$sex))
dat$tumor.subtype = as.numeric(as.factor(dat$tumor.subtype))

# cheat so we have only numerical variables 
# (for the sake of the example):3
X = model.matrix(~0+.,dat) 
p = ncol(X)
n = nrow(X)

# bootstrap selection rates
B = 10
sel.rates = matrix(0,nrow=B,ncol=p)
nms = colnames(X)
colnames(sel.rates) = nms
sel.rates.other.option = sel.rates

for(b in 1:B){
	# prepare bootstrap resample:
	ib = sample(1:n,n,replace=TRUE)
	xb = X[ib,]
	# find redundant variables to remove:
	co = findCorrelation(cor(xb))
	# keep track of variables that are not removed:
	sel.rates[b,-co] = 1
	# Note: alternatively we can use the variable names 
	# to select columns (but then we must use the names of 
	# the variables we want to "keep"):
	keepers = nms[-co]
	sel.rates.other.option[b,keepers] = 1
	# (The above trick is handy if it's easier to extract 
	# variable names rather than variables indices from the 
	# filter output.)
}
identical(sel.rates.other.option,sel.rates)

# In order to recombine selection rates for categorical variables,
# we'd need to automate the following approach:
sel.rates = as.data.frame(sel.rates)
sel.rates$League = pmin(sel.rates[,"LeagueA"]+sel.rates[,"LeagueN"],1)
sel.rates$LeagueA = NULL
sel.rates$LeagueN = NULL
# or:
# sel.rates$LeagueA = pmin(sel.rates[,"LeagueA"]+sel.rates[,"LeagueN"],1)
# sel.rates$LeagueN = NULL

# calculate selection rates:
apply(sel.rates,2,mean) 

# Final word: there may be better approaches, if you find one please share it!
