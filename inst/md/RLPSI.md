## Restrictive Linear Phenotypic Selection Index (RLPSI-KN).

Kempthorne and Nordskog (1959) proposed an index that incorporates restrictions based on a predeteined level and improvement developed the theory of restricted selection index (ISR), which maximizes the genetic progress of some average Phenotypic variables while the other half is maintained without exchange. According to Cunningham {et al. } (1970), it is considered that an IS is restricted if the average allows one or more variables phenotype remains unchanged (no genetic gain) while in other genetic gain is maximized average. In general, the Kempthorne and Nordskog ISR (1959) has similar properties to IS Smith (1936). 

Consider the problem of increasing grain yield in corn while other phenotypic variables remain without exchange. In this situation Kempthorne and Nordskog (1959) (KN) proposed maximize $\rho_{YZ}^2$  incorporating constraints matrix genotypic variance-covariance. Bulmer (1980) showed that Kempthorne and Nordskog results (1959) also obtained maximize ${\theta}'G \beta$. In the derivation of the ISR Kempthorne and Nordskog (1959) follows the approach of Bulmer (1980). 

Suppose we have $q$ variables and desired phenotypic improve the average of $qr$ of them while $r$ of them remain without average exchange, the index $Y=\beta'p$ must maximize  $\rho_{YZ}^2$ under such conditions. The restrictions are introduced by through the matrix $W$ of 0s and 1s, where 1s indicate that mean values of the variables remain unchanged, while the 0s indicate average changes in phenotypic variables, such so that $W'G \beta =0$. Let $C =GW$ , since $\rho_{YZ}$ is invariant under changes of scale, can assume that $\beta 'P \beta = 1$, so $\theta 'G \beta$ should be maximized under the constraint $\beta'C = 0$ and $\beta' P \beta = 1$. 

Let $\tau$ and  $ v'=[v_1 ... v_r]$  Lagrange multipliers. 

Should be maximized 
$\Psi =\theta'G \beta - \tau  (\beta'P \beta-1) -v 'C' \beta$

In deriving $\Psi$ over $\beta$ and equate the result to null vector is $\theta'G-2 \tau P \beta -v'C' =0$

Thus, the vector of coefficients of KN that maximizes $\rho_{YZ}^2$ is  $\beta_{KN}=[I -P^{-1}C (C'P^{-1}C)^{-1}C']\beta_S$. If $A=[ I - P^{-1} C(C'P ^{-1}C)^{-1}C']$, then $\beta_{KN} =A \beta_S$.

```R

library(Rindsel)
datos<-data.frame(read.csv("https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/C1_PSI_05_Phen.csv",header=T,na.strings=c(NA,".","-"))) #Raw data to analized.
datos$REP=as.factor(datos$REP)                                                                    #Transform variables to factor.
datos$Block=as.factor(datos$Block)                                                                #Transform variables to factor.
datos$ENTRY=as.factor(datos$ENTRY)                                                                #Transform variables to factor.
file.wgt<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/weights_C1_PSI.csv"    #name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                         #Selection intensity.
design<-"lattice"                                                                                 #Experimental design.
corr<-FALSE                                                                                       #You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                     #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
one.env<-TRUE                                                                                     #Use FALSE for multienvironment trials.
block.ex<-FALSE                                                                                   #Use FALSE always.
softR<-""                                                                                         #Use "" always.
file.covG<-""                                                                                     #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.

KNIndex(datos,file.wgt,selval,design,corr,out="outextKN.txt",outcsv="outKN.csv",rawdata,one.env,block.ex,softR,file.covG)
file.show("outextKN.txt")
```
[Return to examples](https://github.com/RAngelaPG/RIndSel-R/blob/master/Readme.md)
