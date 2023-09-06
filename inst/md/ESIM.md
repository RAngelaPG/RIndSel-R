## Eigen Selection Index Method (ESIM).

The eigen selection index method (or ESIM for its acronym in English) is a selection method based on the IS Smith (1936) in which uses the theory of singular value decomposition to estimate the vector $\beta$ that maximizes the correlation $\rho_{ZY}$.

In this case the first eigenvector elements ($\beta$) of ${P}^{-1}{G}$ determine the  proportion to each character contributes to IS $Y ={\beta}'p$, also the first eigenvalue of ${P} ^{-1}{G }$ is used in the response to selection, such IS was proposed by Ceron-Rojas {et al}. (2008a, 2008b). 

ESIM theory is directly related to the theory of the canonical correlations in which, according to Muirhead (2005), when a set of variables may be divided so natural two subsets of random variables, and what you want is maximize the correlation between the two subsets of variables, canonical correlation analysis is very efficient. Consider the set of phenotypic variables ($p$) and the set of genotypic variables ($g$) in the context the improvement of plant and animal, in this case must be, indeed, the vector of phenotypic variables, $p$ , and the vector of genotypic variables, $g$ , belong, naturally, two sets of variables. Analysis canonical correlations reduces the correlation between $p$ and $g$ to its simplest form by linear transformations of $p$ and $g$ , ie $Y=\beta'p$ and $Z={\theta'g}$. 

Using a concept similar to that of Kempthorne and Nordskog (1959), Ceron-Rojas {et al}. (2008a, 2008b) maximized the response to selection (Equation 2) to maximize $\rho_{YZ}^2$ . Note that the variances $Y=\beta'p$ and $Z={\theta'g}$ are constant each selection cycle, therefore, the selection of genotypes can be done using $Y=\beta'p$ or $Y/\sqrt{{\beta'}{P\beta}}$. Because of this, by maximizing $\rho_{YZ}^2$ is possible to introduce the constraints $\beta'P \beta=1$ and $\theta'G \theta=1$ so that in ESIM is necessary to maximize $\Phi = (\theta'G \beta)^2 - \mu  (\beta'P \beta -1 ) - \omega  (\theta'G \theta -1 )$ with respect to $\beta$ , $\theta$ , $\mu$ , and $\omega$, where $\beta$ is the vector of coefficients $Y=\beta'p$ , $\theta$ is the vector coefficients of $Z =\theta'g$ , and $\mu$ and $\omega$ are multipliers Lagrange. In ESIM the values of $\theta$ does not necessarily economic weights. 

When $\Phi$ is derived with respect to $\beta$ and $\theta$ And the result is equal to zero vector, we have $(\theta 'G \beta)^{2}$, $G\theta - \mu P \beta = 0$, $(\theta 'G \beta)^{2}$, $G \beta - \omega G \theta =0$.
 
By the constraints $\beta'P \beta=1$ and $\theta'G \theta=1$ ,when the last equation  is multiplied by $\beta'$ and the next equation is multiplied by $\theta'$ ,the result is $({\theta}'G \beta)^2 = \omega= \mu$.Therefore, $\mu$ maximizes $\rho_{YZ}^2$ under the constraints $\beta' P \beta =1$ and $\theta' G \theta=1$.

The next problem is to determine the vector $\beta$ that allows building the IS $Y =\beta'p$ that has maximum correlation with $Z=\theta'g$. According to Ceron-Rojas {et al}. (2008a, 2008b), the $\beta$ required is the solution to the following equation 

$( Q - \mu I)\beta = 0$

where $Q= P ^{-1} G$ . Thus, in ESIM , the value that maximizes $\rho_{YZ}^2$ under the constraints $\beta'P \beta=1$ and $\theta'G \theta=1$ is the first eigenvalue ($\mu$) matrix $Q$ , and the vector allows to construct $Y=\beta'p$ (with maximum correlation with $Z=\theta'g$) is the first eigenvector ($\beta$) Matrix $Q$. 

```R

library(Rindsel)
datos<-read.csv("https://github.com/RAngelaPG/RIndSel-R/blob/master/data/C1_PSI_05_Phen.csv",header=T,na.strings=c(NA,"."."-")) #Raw data to analized.
file.wgt<-"https://github.com/RAngelaPG/RIndSel-R/blob/master/data/weigth_C1_PSI.csv")       #name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                    #Selection intensity.
design<-"lattice"                                                                            #Experimental design.
corr<-FALSE                                                                                  #You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
one.env<-TRUE                                                                                #Use FALSE for multienvironment trials.
block.ex<-FALSE                                                                              #Use FALSE always.
softR<-""                                                                                    #Use "" always.
file.covG<-""                                                                                #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.

ESIMIndex(datos,file.wgt,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",rawdata,one.env,block.ex,softR,file.covG)

```
[Return to examples](https://github.com/RAngelaPG/RIndSel-R/blob/master/Readme.md)


