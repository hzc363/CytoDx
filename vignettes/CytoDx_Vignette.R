## ---- out.width = "500px",echo=FALSE-------------------------------------
knitr::include_graphics("tssm_intro.jpg")

## ---- eval=FALSE---------------------------------------------------------
#  BiocInstaller::biocLite("CytoDx")

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("hzc363/CytoDx")

## ---- results='asis',message=FALSE---------------------------------------
library(CytoDx)

# Find data in CytoDx package
path <- system.file("extdata",package="CytoDx")

# read the ground truth
fcs_info  <- read.csv(file.path(path,"fcs_info.csv"))

# print out the ground truth
knitr::kable(fcs_info)

## ---- results='asis',message=FALSE---------------------------------------
# Find the training data
train_info <- subset(fcs_info,fcs_info$dataset=="train")

# Specify the path to the cytometry files
fn <- file.path(path,train_info$fcsName)

# Read cytometry files using fcs2DF function
train_data <- fcs2DF(fcsFiles=fn,
                    y=train_info$Label,
                    assay="FCM",
                    b=1/150,
                    excludeTransformParameters=
                      c("FSC-A","FSC-W","FSC-H","Time"))

## ---- results='asis',message=FALSE---------------------------------------
# Perfroms rank transformation
x_train <- pRank(x=train_data[,1:7],xSample=train_data$xSample)

# Convert data frame into matrix. Here we included the 2-way interactions.
x_train <- model.matrix(~.*.,x_train)

## ---- results='asis',message=FALSE,warning=FALSE-------------------------
# Build predictive model using the CytoDx.fit function
fit <- CytoDx.fit(x=x_train,
               y=(train_data$y=="aml"),
               xSample=train_data$xSample,
               family = "binomial",
               reg = FALSE)

## ---- results='asis',message=FALSE---------------------------------------
# Find testing data
test_info <- subset(fcs_info,fcs_info$dataset=="test")

# Specify the path to cytometry files
fn <- file.path(path,test_info$fcsName)

# Read cytometry files using fcs2DF function
test_data <- fcs2DF(fcsFiles=fn,
                    y=NULL,
                    assay="FCM",
                    b=1/150,
                    excludeTransformParameters=
                      c("FSC-A","FSC-W","FSC-H","Time"))
# Perfroms rank transformation
x_test <- pRank(x=test_data[,1:7],xSample=test_data$xSample)

# Convert data frame into matrix. Here we included the 2-way interactions.
x_test <- model.matrix(~.*.,x_test)

## ---- results='asis',message=FALSE---------------------------------------
# Predict AML using CytoDx.ped function
pred <- CytoDx.pred(fit,xNew=x_test,xSampleNew=test_data$xSample)

## ---- results='asis',message=FALSE,fig.width = 5-------------------------
# Cmbine prediction and truth
result <- data.frame("Truth"=test_info$Label,
                    "Prob"=pred$xNew.Pred.sample$y.Pred.s0)

# Plot the prediction
stripchart(result$Prob~result$Truth, jitter = 0.1,
           vertical = TRUE, method = "jitter", pch = 20,
           xlab="Truth",ylab="Predicted Prob of AML")

## ---- results='asis',message=FALSE,fig.width = 5-------------------------
# Use decision tree to find the cell subsets that are associated the AML.
TG <- treeGate(P = fit$train.Data.cell$y.Pred.s0,
              x= train_data[,1:7])


## ---- results='asis',message=FALSE,fig.width = 5-------------------------
sessionInfo()


