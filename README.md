## Introduction

CytoDx is a method that predicts clinical outcomes using single cell data without the need of cell gating. It first predicts the association between each cell and the outcome using a linear statistical model. The cell level predictions are then averaged within each sample to represent the sample level predictor. A second sample level model is used to make prediction at the sample level.  


## Installation

To install the MetaCyto package, please run the following code:
```
library("devtools")
install_github("hzc363/CytoDx")
```

## Example
Here we provide a simple simulated example. Please see the vignette for an example that applies tssm to diagnos accute myeloid lymphoma using flow cytometry data. 

```
library(tssm)

# simulate 10 samples of class A 
A = data.frame("marker1"=c(rnorm(1000),rnorm(1000,mean=10)),
                "marker2"=c(rnorm(1000),rnorm(1000,mean=10)),
                "sample"=sample(1:10,2000,replace = T),
                "y"="A")

# simulate 10 samples of class B
B = data.frame("marker1"=c(rnorm(1500),rnorm(500,mean=10)),
               "marker2"=c(rnorm(1500),rnorm(500,mean=10)),
               "sample"=sample(11:20,2000,replace = T),
               "y"="B")             

# build CytoDx model
dat = rbind(A,B)
fit = CytoDx.fit(x = as.matrix(dat[,1:2]),
               y = (dat$y=="A"),
               xSample = dat$sample,
               family = "binomial",
               reg = F)

# predict 
pred = CytoDx.pred(fit = fit,
                 xNew = as.matrix(dat[,1:2]),
                 xSampleNew = dat$sample)

# plot probability
result = data.frame("Truth"=rep(c("A","B"),each=10),
                    "Prob"=pred$xNew.Pred.sample$y.Pred.s0)
stripchart(result$Prob~result$Truth, jitter = 0.1,
           vertical = TRUE, method = "jitter", pch = 20,
           xlab="Truth",ylab="Predicted Prob of A")
```
