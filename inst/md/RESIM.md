## Restrictive Eigen Selection Index Method (RESIM).

RESIM is a generalization of ESIM to the case in which are incorporated restrictions on $Y ={\beta'p}$ , the procedure for set the $\beta$ is, therefore, very similar to ESIM and is based on IS restricted Kempthorne and Nordskog (1959). So, suppose you have $q$ variables and phenotypic want to improve the means of $q-r$ of them while $r$ of them average remain unchanged, the index $Y=\beta'p$ should maximize $\rho_{YZ}^2$ under conditions tale. Again, the restrictions are introduced by the matrix $W$ of 0s and 1s, where 1s indicate that the mean values of the variables Phenotypic remain unchanged, while indicating 0s average changes in phenotypic variables, such that $W'G \beta =0$. Let $C=GW$, since $\rho_{YZ}$ is invariant under changes of scale can be assumed which $\beta'P \beta=1$ and $\theta'G \theta=1$ , so $\rho_{YZ}^2$ be maximized under the constraints $\beta'C=0$, $\beta'P \beta=1$ and $\theta'G \theta=1$. Let $\mu$, $\omega$ and $v'=[v_1 ... v_r]$ Lagrange multipliers. It should maximized $\Phi=(\theta'G\beta)^2 - \mu (\beta'P\beta-1) - \omega (\theta'G\theta-1)-v'C'\beta$

In deriving $\Phi$ over $\beta$, $\theta$, $\mu$, $\omega$, and $v'=[v_1 ... v_r]$, Ceron-Rojas {et al}. (2008a) showed that $(Q-\omega I)\beta=0$

where now $Q=[I-P^{-1}C(C'P^{-1}C)^{-1}C']P^{-1}G$. Thus, RESIM, the value that maximizes $\rho_{YZ}^2$ under the constraint $\beta'P\beta=1$ , $\theta'G\theta=1$ and $\beta'C=0$ is the first eigenvalue ($\omega$) of the matrix $Q$ , and the vector lets build resim $Y=\beta'p$ (with maximum correlation with $Z=\theta'g$) is the first eigenvector ($\beta$) matrix $Q$.

```R

library(Rindsel)
datos<-read.csv("https://github.com/RAngelaPG/RIndSel-R/blob/master/data/C1_PSI_05_Phen.csv",header=T,na.strings=c(NA,"."."-")) #Raw data to analized.
file.wgt<-"https://github.com/RAngelaPG/RIndSel-R/blob/master/data/weigth_C1_PSI.csv")   #name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                    #Selection intensity.
design<-"lattice"                                                                            #Experimental design.
corr<-FALSE                                                                                  #You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
one.env<-TRUE                                                                                #Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                              #Use FALSE always.
softR<-""                                                                                    #Use "" always.
file.covG<-""                                                                                #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.

RESIMIndex(datos,file.wgt,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",rawdata,one.env,block.ex,softR,file.covG)

```
[Return to examples](https://github.com/RAngelaPG/RIndSel-R/blob/master/Readme.md)
