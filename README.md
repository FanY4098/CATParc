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
```
After obtaining the test statistics (i.e., coupling scores) for each pair of positions, the AUC can be computed using the corresponding distance file provided in the data/ folder.

# Implementation of Competing Methods

This note is intended to facilitate reproducibility of the results reported in Section 5.1 of our manuscript, “Model-Free Inference for Characterizing Protein Mutations through a Coevolutionary Lens.” In this section, we compare our method with several existing approaches, including PSICOV, GREMLIN, plmDCA, mfDCA, and GaussDCA, across six protein families. Note that PSICOV and GREMLIN use input datasets in .txt format, whereas plmDCA, mfDCA, and GaussDCA, take datasets in .fasta format. All datasets are available in the data/ folder on GitHub. In addition, for a fair comparison with our method, we run these methods without sequence reweighting or APC correction.

For illustration purposes, we use PF01037 as an example. The working directory is assumed to be ./ and should be adjusted accordingly when running the commands. Output files are saved with names ending in result.txt or result.csv.

PSICOV

PSICOV can be installed from the official GitHub repository at https://github.com/psipred/psicov. After installation, run the following command in the terminal to obtain the scores for each pair of positions:

```
./psicov -d 0.03 -j 1 -g 1 -l PF01037.txt > PF01037.psicov.result.txt
```

GREMLIN

GREMLIN can be installed from the official GitHub repository at https://github.com/sokrypton/GREMLIN/tree/master. After installation, run the following command in MATLAB to obtain the scores for each pair of positions.

```
gremlin(‘./PF01037.txt’, ‘./PF01037.gremlin.result.txt’, 'apc', 0, ‘reweight’, 0)
```

plmDCA

plmDCA can be installed from the GitHub repository at https://github.com/debbiemarkslab/plmc. The current pipeline does not directly output coupling scores without APC correction. Therefore, we first save the estimated model parameters and then manually compute the coupling scores for each pair of positions.

After installation, run the following command in the terminal to save the estimated parameters to a .params file. Note that, to run plmDCA without sequence reweighting, we set -s to the total number of sequences in the MSA (e.g., 35725 sequences for PF01037).

```
./plmc -o ./PF01037.params -t 1 -s 35725 -m 100 ./PF01037.fasta
```

Then run the following commands in MATLAB to obtain the final coupling scores (saved as results.txt).

```
cd(‘./plmc-master/scripts’)

params=read_params(‘./PF01037.params’)

L = size(params.Jij, 1);    % number of positions
q = size(params.Jij, 3);    % number of amino acid states

results = [];  % to store [i, j, score]

for i = 1:L
    for j = i+1:L
        Jij = squeeze(params.Jij(i, j, :, :));  % q x q coupling matrix
        score = norm(Jij, 'fro');               % Frobenius norm
        results(end+1, :) = [i, j, score];      % append to result
    end
end

writematrix(results, ‘./PF01037.plmdca.result.txt', 'Delimiter', 'tab')
```

mfDCA

mfDCA can be installed from the GitHub repository at https://github.com/utdal/py-mfdca. After installation, run the following Python commands to obtain the desired results.

```
from dca.dca_class import dca
import numpy as np
import time



pf = dca('/Users/fan/Documents/PITT/#Research_Pitt/Protein Project/plmc-master/bin/PF01037.fasta’)
pf.mean_field(theta=1)  # disable reweighting

L = pf.couplings.shape[0]         # number of positions
scores = []
for i in range(L):
    for j in range(i+1, L):
        Jij = pf.couplings[i, j, :, :]      # q × q coupling matrix
        score = np.linalg.norm(Jij, 'fro')   # Frobenius norm
        scores.append((i+1, j+1, score))     # 1-based indexing if desired
np.savetxt(‘./PF01037.mfdca.result.csv’,scores, delimiter=",")
```

GaussDCA

GaussDCA can be installed from the GitHub repository at https://github.com/carlobaldassi/GaussDCA.matlab. After installation, run the following MATLAB commands to obtain the desired results.

```
R= gDCA(‘./PF01037.fasta', ’max_gap_fraction’, 1, ‘score’, ‘frob’, 'theta', 0, ’min_separation’, 1, ‘apply_apc’, false)
T = array2table(R, 'VariableNames', {'i', 'j', 'score'})

writetable(T, ‘./PF01037.gaussdca.result.csv’)
```

