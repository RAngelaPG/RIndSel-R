## Molecular Eigen Selection Index Method (MESIM).

MESIM is a generalization of ESIM to the case in which is incorporated MM information to IS, similarly as in the IS Lande and Thompson (1990). Following the basic idea of Kempthorne and Nordskog (1959), Ceron-Rojas {et al}. (2008b), maximized the correlation between $Y_M=\beta'_p p + \beta'_s s=\beta'_M p_ps$ and $Z_M=\theta'_1g+\theta'_2s=\theta'_M g{gs}$ , $\rho{YM ZM}^2$, where $s$ is the vector of records of the additive effects of QTLs associated with MM, $p'{ps}=[p's']$, $g'{gs}=[g's']$, $\beta'_M=[\beta'_P\beta'_s]$, and $\theta'_M=[\theta'_1\theta'_2]$.

In the illustrative example of Lande and Thompson (1990), $Y_M=\beta_p*p+\beta_s*s$, and $Z_{M}=\theta_{p}*G_{P}+\theta_{s}*s$ , where $p$ denotes the plant height and $s=x_1 \alpha_1 + x_2 \alpha_2 + ... + x_5 \alpha_5$. 

Again, because $\rho_{Y_{M}Z_{M}}^2$  is invariant to changes in is maximized scale $\rho_{Y_M Z_M ^2}$ under the constraints $\beta'_MT\beta_M=1$ and $\beta'_{M}*T*\beta_{M}=1$ and $\theta'_{M}K\theta_{M}=1$ and $\theta'_M*K*\theta_M=1$, so in MESIM, it is necessary maximize 

$$\Phi=(\theta'_MK\beta_M)^2 - \mu (\beta'_MT\beta_M-1) - \omega(\theta'_MK\theta_M-1)$$

With respect to $\beta_M$, $\theta_M$, $\mu$ , and $\omega$, where $\beta_M$ is the vector of coefficients MESIM, $\theta_M$ is the vector of coefficients $Z_{M}=\theta'_{M}g_{gs}$ and $\mu$ and $\omega$ are Lagrange multipliers. Ceron-Rojas {et al}. (2008b) found that the solution is equal 

$$(Q-\mu I)\beta_M=0$$

where $Q=T^{-1}K$. Thus, in MESIM, the value that maximizes $\rho_{Y_M Z_M}^2$ under the constraints $\beta'_{M}T\beta_{M}=1$ and $\theta'_{M}K\theta_{M}=1$ is the first eigenvalue ($\mu$) Matrix $Q$, and the vector for building $Y_M$ (With maximum correlation with $Z_{M}=\theta'_{M}g_{gs})$ is the first eigenvector ($\beta_M$) matrix $Q$.

```R

library(Rindsel)
datos<-read.csv("https://github.com/RAngelaPG/RIndSel-R/blob/master/data/C1_PSI_05_Phen.csv",header=T,na.strings=c(NA,"."."-")) #Raw data to analized.
file.wgt<-"https://github.com/RAngelaPG/RIndSel-R/blob/master/data/weigth_C1_PSI.csv")             #name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                          #Selection intensity.
design<-"lattice"                                                                                  #Experimental design.
corr<-FALSE                                                                                        #You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                      #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
file_nameMARK<-"https://github.com/RAngelaPG/RIndSel-R/blob/master/data/C1_PSI_S2_05_Haplo.csv")   #name of the file markers information.
one.env<-TRUE                                                                                      #Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                                    #Use FALSE always.
softR<-""                                                                                          #Use "" always.
file.covG<-""                                                                                      #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.

MESIMIndex(datos,file.wgt,selval,design,corr,out="outextLT.txt",outcsv="outLT.csv",rawdata,file_nameMARK,one.env,block.ex,softR,file.covG)

```
[Return to examples](https://github.com/RAngelaPG/RIndSel-R/blob/master/Readme.md)
