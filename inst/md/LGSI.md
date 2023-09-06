## Linear Genomic Selection Index (LGSI).

In a similar manner as the linear phenotypic selection index (LPSI), the objective of the LGSI is to predict the net genetic merit $H=w'g$, where $g'=(g_1,g_2,...,g_t)$ ($t$ number of traits) is a vector of unobservable true breeding values and $w'=(w_1,w_2,...,w_t)$ is a vector of economic weights. Suppose that the genomic breeding values $\gamma_i=Xu_i$ are known; then, the LGSI can be written as
```math
\tag{1}
\begin{equation}I_G=\beta' \gamma \end{equation}
```
where $\beta$ is an unknown vector of weights.

The main advantage of the LGSI over the LPSI lies in the possibility of reducing the intervals between selection cycles ($L_G$) by more than two thirds (Lorenz et al. 2011); thus, this parameter should be incorporated into the LGSI selection response and the expected genetic gain per trait to reflect the main advantage of the LGSI over the LPSI and the other indices. In the LPSI context we wrote the selection response as $R_I=k_I\sigma_H \rho_{HI}$, while the LGSI selection response can be written as
```math
\tag{2}
\begin{equation}R_{I_G}=\frac{k_I \sigma_{HI_G}}{L_G \sigma^2_{I_G}}=\frac{k_I}{L_G}\sigma_H \rho_{HI_G}\end{equation}
```
where $k_I$ is the standardized selection differential (or selection intensity) associated with LGSI, $\sigma_{HI_G}$ is the covariance between $H=w'g$ and LGSI, $\sigma^2_{I_G}$ is the variance of LGSI, $\sigma_H$ is the standard deviation of $H$, $\rho_{HI_G}$ is the correlation between $H$ and LGSI, and $L_G$ denotes the intervals between selection cycles. 
	
Let $C$ and $\Gamma$ be matrices of covariance of the breeding values ($g$) and of the genomic breeding values ($\gamma$), respectively; then, the correlation between $H=w'g$ and $I_G=\beta' \gamma$ can be written as
```math
\tag{3}
\begin{equation}\rho_{HI_G}=\frac{w' \Gamma \beta}{\sqrt{w'Cw}\sqrt{\beta' \Gamma \beta}}\end{equation}
```
where $w' \Gamma \beta= \sigma_{HI_G}$ is the covariance between $H=w'g$ and $I_G=\beta' \gamma$ , $\sigma_H=\sqrt{w'Cw}$ is the standard deviation of the variance of $H=w'g$, and $\sigma_{I_G}=\sqrt{\beta' \Gamma \beta}$ is the standard deviation of the variance of $I_G=\beta' \gamma$.

###The maximize LGSI parameters.

To maximize the genomic selection response (Equation 2), suppose that $k_I$, $\sigma_H$ and $L_G$ are fixed and take the derivative of the natural logarithm (ln) of the correlation between $H$ and $I_G)$ (Equation 3) with respect to vector $\beta$, equate the result of the derivative to the null vector, and isolate $\beta$, i.e.,
```math
\tag{4}
\begin{equation}\frac{d}{d\beta}ln \rho_{HI_G}=\frac{d}{d\beta}ln \left(\frac{w'\Gamma \beta}{\sqrt{w'Cw}\sqrt{\beta' \Gamma \beta}} \right)=0\end{equation}
```
The result is $\beta=sw$, where $s=\frac{\beta' \Gamma \beta}{w' \Gamma \beta}$ is a proportional constant that does not affect the maximum value of $\rho_{HI_G}$ because this is invariant to the scale change; then assuming that $\beta=w$, the maximized LGSI selection response can be written as
```math
\tag{5}
\begin{equation}R_{I_G}=\frac{k_I}{L_G}\sqrt{w' \Gamma w}\end{equation}
```
Hereafter we will refer to the LGSI genomic selection response as that of Equation (5). 

Also, because $\beta=w$, Equation (3) can be written as
```math
\tag{6}
\begin{equation}\rho_{HI_G}=\frac{\sqrt{w' \Gamma w}}{\sqrt{w' C w}}=\frac{\sigma_{I_G}}{\sigma_H}\end{equation}
```
which is the maximized correlation between $H=w'g$ and $I_G=\beta' \gamma$, or LGSI accuracy; $\sigma_H=\sqrt{w'Cw}$ is the standard deviation of the variance of $H$, and $\sigma_{I_G}=\sqrt{\beta' \Gamma \beta}$ is the standard deviation of the variance of $I_G$.

The LGSI expected genetic gain per trait ($E_{I_G}$) can be written as
```math
\tag{7}
\begin{equation}E_{I_G}=\frac{k_I \Gamma w}{L_G \sqrt{w' \Gamma w}}\end{equation}
```
All the terms in Equation (7) were previously defined.

Let $\lambda_G=\frac{\rho_{HI_G}}{\rho_{HI}}$ be LGSI efficiency vs. LPSI efficiency to predict the net genetic merit, where $\rho_{HI_G}$ is the LGSI accuracy and $\rho_{HI}$ the LPSI accuracy; in percentage terms, LGSI efficiency vs. LPSI efficiency for each selection cycle can be written as
```math
\tag{8}
\begin{equation}p_G=100(\lambda_G-1)\end{equation}
```
According to Equation (8), if $p_G>0$, LGSI efficiency will be greater than LPSI efficiency; if $p_G=0$, the efficiency of both selection indices will be equal, and if $p_G<0$, LPSI will be more efficient than LGSI for predicting $H=w'g$.

Equation (8) is useful for measuring LGSI efficiency in terms of accuracy when predicting the net genetic merit ($H=w'g$), while the Technow et al. (2013) inequality measures LGSI efficiency in terms of the time needed to complete one selection cycle. In the context of LGSI and LPSI, the Technow inequality can be written as
```math
\tag{9}
\begin{equation}L_G<\frac{\rho_{HI_G}}{h_I}L_P\end{equation}
```
where $L_G$ and $L_P$ denote the time required to complete one selection cycle for LGSI and LPSI, respectively; $\rho_{HI_G}$ is the LGSI accuracy, and $h_I$ is the square root of the heritability (Lin and Allaire 1977; Nordskog 1978) of LPSI, which can be denoted as $h_I=\sqrt{\frac{b'Cb}{b'Pb}}$. Then, assuming that the selection intensity is the same for both selection indices, if Equation (9) is true, LGSI will be more efficient than LPSI per unit of time.

###Statistical LGSI properties.
Assuming that $H$ and $I_G$ have joint bivariate normal distribution and that $\Gamma$, $C$ and $w$ are known, the LGSI has the following properties:

* The variance of $I_G$ ($\sigma^2_{I_G}$) and the covariance between $H$ and $I_G$ ($\sigma_{HI_G}$) are equal, i.e., $\sigma^2_{I_G}=\sigma_{HI_G}$.
* The maximized correlation between $H$ and $I_G$ (or LGSI accuracy) is equal to $\rho_{HI_G}=\frac{\sigma_{I_G}}{\sigma_H}$, where $\sigma_{I_G}$ is the standard deviation of $\sigma^2_{I_G}$ and $\sigma_H$ is the standard deviation of the variance of $H$ ($\sigma^2_H$).
* The variance of the predicted error, $Var(H-I_G)=(1-\rho^2_{HI_G})\sigma^2_H$, is minimal. Note that $Var(H-I_G)=\sigma^2_{I_G}+\sigma^2_H-2\sigma_{HI_G}$, and when $\beta=w$, $\sigma^2_{I_G}=\sigma_{HI_G}$, from where $Var(H-I_G)=\sigma^2_H-\sigma^2_{I_G}=(1-\rho^2_{HI_G})\sigma^2_H$ is minimal.
* The total variance of $H$ explained by $I_G$ is $\sigma^2_{I_G}=\rho^2_{HI_G}\sigma^2_H$. It is evident that if $\rho_{HI_G}=1$, $\sigma^2_{I_G}=\sigma^2_H$, and if $\rho_{HI_G}=0$, $\sigma^2_{I_G}=0$. That is, the variance of $H$ explained by $I_G$ is proportional to $\rho_{HI_G}$, and when $\rho_{HI_G}$ is close to 1, $\sigma^2_{I_G}$ is close to $\sigma^2_{H}$; if $\rho_{HI_G}$ is close to 0, $\sigma^2_{I_G}$   is close to 0.

The LGSI properties described in points 1 to 4 of this subsection are the same as the LPSI properties described. This corroborates that the LGSI is an application of the LPSI theory to the GS context.

```R
library(Rindsel)
datos<-read.csv("https://github.com/RAngelaPG/RIndSel-R/blob/master/data/Data_Phenotypes_LGSI.csv",header=T,na.strings=c(NA,"."."-")) #Raw data to analized.
file.wgt<-"https://github.com/RAngelaPG/RIndSel-R/blob/master/data/weigths_LGSI.csv")   		#name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                    		#Selection intensity.
design<-"lattice"                                                                            		#Experimental design.
corr<-FALSE                                                                                  		#You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                		#By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
file_nameMARK<-"https://github.com/RAngelaPG/RIndSel-R/blob/master/data/Training population_LGSI.csv")  #name of the file training markers information.
one.env<-TRUE                                                                                		#Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                              		#Use FALSE always.
softR<-""                                                                                    		#Use "" always.
file.covG<-""                                                                                		#When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.
file_nameTMARK<-"https://github.com/RAngelaPG/RIndSel-R/blob/master/data/Testing population_LGSI.csv")	#name of the file testing markers information.
LG<-1													#Interval between selection cycles.

LGSI(datos,file.wgt,selval,design,corr,method,out="outextLGSI.txt",outcsv="outLGSI.csv",rawdata,file_nameMARK,one.env,block.ex,softR,file.covG,file_nameTMARK,LG)

```
[Return to examples](https://github.com/RAngelaPG/RIndSel-R/blob/master/Readme.md)
