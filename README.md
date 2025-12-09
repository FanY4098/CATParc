# CATParc
The package implements a spectrum-based statistic to test if two positions along the protein sequence are partially correlated or not.

# Installation
To install the `CATParc` package, you will first need to install `devtools` package and then execute the following code: 
```
devtools::install_github('FanY4098/CATParc')
```

# Main Functions
There are four main functions in this package. (1)`OneHot_Protein` transforms the original MSA data into standardized binary vectors; (2) `ModFit` fits regularized regression model based on the transformed data to get the residuals; (3) Based on the fitted residuals, `TestCCA` employs a spectrum based test statistics to infer the partial correlation between two positions of the protein sequences, while `TestLS` uses the L2 and Sup norm based statistics to do the inference. You can always use the following command to see more details:
```
library(CATParc)
?OneHot_Protein
?ModFit
?TestCCA
?TestLS
```

# Application on PF00502
We illustrate the use of our package by testing if positions 1 and 2 of PF00502 form a contact (i.e., to test if the partial correlation is 0 or not.). It can be executed using the following codes.

```
data(SixProteinFamilies)

##data preprocessing
pre=OneHot_Protein(PF00502)
data=pre$data1hot
Gsz=pre$Gsz

##fit models to get residuals
Gnm=length(Gsz)
allend=cumsum(Gsz) #ending position for each group
allstart=allend+1
allstart=c(1,allstart[-Gnm]) #starting position for each group

y1=as.matrix(data[,allstart[1]:allend[1]])
y2=as.matrix(data[,allstart[2]:allend[2]])
x=data[,-c(allstart[1]:allend[1],allstart[2]:allend[2])]
xGsz=Gsz[-c(1,2)]


mod1=ModFit(y=y1,x=x,xGsz=xGsz,lam1=0,lamG=0.07*(sqrt(xGsz*dim(y1)[2]/dim(x)[1])+sqrt(2*log(length(xGsz))/dim(x)[1])))
res1=mod1$res
mod2=ModFit(y=y2,x=x,xGsz=xGsz,lam1=0,lamG=0.07*(sqrt(xGsz*dim(y2)[2]/dim(x)[1])+sqrt(2*log(length(xGsz))/dim(x)[1])))
res2=mod2$res

##get p-values from different tests
CCA=TestCCA(res1,res2)$pval
L2=TestLS(res1,res2)$pval_chisq
SUP=TestLS(res1,res2)$pval_sup
