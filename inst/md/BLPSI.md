## Base Linear Phenotypic Selection Index (BLPSI)

To derive the LPSI theory, we assumed that the phenotypic ($P$) and the genotypic ($G$) covariance matrix, and the vector of economic values ($w$) are known. However, $P$, $G$ and $w$ are generally unknown and it is necessary to estimate them. There are many methods for estimating $P$ and $G$ (Lynch and Walsh 1998) and $w$ (Cotterill and Jackson 1985; Magnussen 1990). However, when the estimator of $P$ ($\hat{P}$) is not positive definite (all eigenvalues positive) or the estimator of $G$ ($\hat{G}$) is not positive semidefinite (no negative eigenvalues), the estimator of $G$ ($\hat{G}$) could be biased. In this case, the base linear phenotypic selection index (BLPSI):

```math
\tag{1}
\begin{equation}
I_{B}=w'y
\end{equation}
```

may be a better predictor of $H=w'g$ than the estimated LPSI $\hat{I}=\hat{b}'y$ (Williams 1962a; Lin 1978) if indeed the vector of economic values   is known. Many authors (Williams 1962b; Harris 1964; Hayes and Hill 1980, 1981) have investigated the influence of parameter estimation errors on LPSI accuracy and concluded that those errors affect the accuracy of $\hat{I}=\hat{b}'y$ when the accuracy of $\hat{P}$ and $\hat{G}$ is low. If vector $w$ values are known, the BLPSI has certain advantages because of its simplicity and its freedom from parameter estimation errors (Lin 1978). Williams (1962a) pointed out that the BLPSI is superior to $\hat{I}=\hat{b}'y$ unless a large amount of data is available for estimating $P$ and $G$.

There are some problems associated with the BLPSI. For example, what is the BLPSI selection response and the BLPSI expected genetic gains per trait when no data are available for estimating $P$ and $G$? The BLPSI is a better selection index than the standard LPSI only if the correlation between the BLPSI and the net genetic merit is higher than the correlation between the LPSI and the net genetic merit (Hazel 1943). But if estimations of $P$ and $G$ are not available, how can we obtain the correlation between the base index and the net genetic merit? Williams (1962b) pointed out that the correlation between the BLPSI and $H=w'g$ can be written as

```math
\tag{2}
\begin{equation}
\rho_{HI_{B}}=\sqrt{\frac{w'Gw}{w'Pw}}
\end{equation}
```
and indicated that the ratio $\frac{\rho_{HI_{B}}}{\rho_{HI}}$ can be used to compare LPSI efficiency vs BLPSI efficiency; however,  in the latter case we at least need to know the estimates of $P$ and $G$, i.e., $\hat{P}$ and $\hat{G}$.

In addition, Equation (1) is only an assumption, not a result, and implies that $P$ and $G$ are the same. That is, $b=P^{-1}Gw=w$ only when $P=G$, which indicates that the BLPSI is a special case of the LPSI. Thus to obtain the selection response and the expected genetic gains per trait of the BLPSI, we need some information about $P$ and $G$. Assuming that BLPSI is indeed a particular case of the LPSI, the BLPSI selection response and the BLPSI expected genetic gains per trait could be written as 

```math
\tag{3}
\begin{equation}
R_B=k_I \sqrt{w'Pw}
\end{equation}
```
and 
```math
\tag{4}
\begin{equation}
E_B=k_I \frac{Gw}{\sqrt{w'Pw}}
\end{equation}
```
respectively. The parameters of Equations (3) and (4) were defined earlier.

There are additional implications if $b=P^{-1}Gw$. For example, if $P=G$, then $\rho_{HI_B}=\sqrt{\frac{w'Gw}{w'Pw}}$ and BLPSI heritability $h^2_{I_B}=\frac{w'Gw}{w'Pw}$ are equal to 1. However, in practice, the estimated values $\hat{\rho_{HI_B}}$ are usually lower than the estimated values $\hat{\rho_{HI}}$.

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

BLPSI(datos,file.wgt,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",rawdata,one.env,block.ex,softR,file.covG)
file.show("outextBLPSI.txt")
```

[Return to examples](https://github.com/RAngelaPG/RIndSel-R/blob/master/Readme.md)
