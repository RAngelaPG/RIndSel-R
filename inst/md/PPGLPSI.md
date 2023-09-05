#### Predeterminated Proportional Gain Linear Phenotypic Selection Index (PPG-LPSI).

This index is called predetermined proportional gain linear phenotypic selection index (PPG-LPSI) because the breeder pre-sets optimal levels for certain traits before the selection is carried out. Some of the main objectives of the PPG-LPSI are to optimize the expected genetic gain per trait, predict the net genetic merit $H=w'g$ and select the individuals with the highest net genetic merit values as parents of the next generation. The PPG-LPSI allows imposing restrictions different from zero on the expected genetic gains of some traits, while other traits increase (or decrease) their expected genetic gains without imposing any restrictions. The PPG-LPSI solves the LPSI equations subject to the condition that the covariance between the LPSI and some linear functions of the genotypes involved be equal to a vector of predetermined constants or genetic gains defined by the breeder.

Let $d'=(d_1, d_2,...,d_r)$ be a vector $rx1$ of the predetermined proportional gains and assume that $\mu_q$ is the population mean of the $q^{th}$ trait before selection. One objective could be to change $\mu_q$ to $\mu_q+d_q$, where $d_q$ is a predetermined change in $\mu_q$ (in the RLPSI, $d_q=0$, $q=1,2,...,r$, where $r$ is the number of predetermined proportional gains). We can solve this problem in a similar manner as we did in the RLPSI. That is, minimizing the mean squared difference between $I$ and $H$ under the restriction $D'U'Gb=0$, where 

$D'=\left( \begin{array}{ccccc}
d_r & 0 & \cdots & 0 & -d_1 \\
0 & d_r & \cdots & 0 & -d_2 \\
\vdots & \vdots & \ddots & \vdots & \vdots\\
0 & 0 & \cdots & d_r & -d_{r-1} \end{array} \right)$ 

is a Mallard (1972) matrix $(r-1)xr$ of predetermined proportional gains, $d_q$ (q=1,2,...,r) is the $q^{th}$ element of vector $d'$, $U'$ is the RLPSI matrix of restrictions of 1’s and 0’s, $G$ is the covariance matrix of genotypic values and $b$ is the LPSI vector of coefficients. Also, it is possible to minimize $E[(H-1)^2]$ under the restriction $U'Gb=\theta d$ (Tallis 1985), where $\theta$ is a proportionality constant. Both approaches are very similar but the equations obtained when introducing the $D'U'Gb=0$ restriction are simpler than when introducing $U'Gb=\theta d$ restrictions into the process of minimizing $E[(H-1)^2]$.

### The maximized PPG-LPSI parameters}
Let $M'=D'C'$ be the Mallard (1972) matrix of predetermined restrictions, where $C'=U'G$. Under the restriction $M'b=0$, we can minimize $E[(H-1)^2]$ assuming that $P$, $G$, $U'$, $D'$ and $w$ are known; that is, we need to minimize the function $\Phi(b,v)=b'Pb+w'Gw-2w'Gb+2v'M'b$ with respect to vectors $b$ and $v'=(v_1, v_2,...,v_{r-1})$, where $v$ is a vector of Lagrange multipliers. Equation (1) derivative results from $b$ and $v$, i.e.,

\[\left( \begin{array}{cc}
P & M \\
M' & 0 
\end{array} \right)
\left( \begin{array}{c}
b \\
v
\end{array} \right)=
\left( \begin{array}{c}
Gw\\
0 
\end{array} \right)\]

from where the vector that minimizes $E[(H-1)^2]$ under the restriction $M'b=0$ is $b_M=K_Mb$ where $K_M=[I_t-Q_M]$, $Q_M=P^{-1}M(M'P^{-1}M)^{-1}M'=P^{-1}CD(D'C'P^{-1}CD)^{-1}D'C'$ and $I_t$ is an identity matrix of size $txt$. When $D=U$, $b_M=b_R$ (the RLPSI vector of coefficients), and when $D=U$ and $U'$ is a null matrix, $b_M=b$ (the LPSI vector of coefficients). Thus, the Mallard (1972) index is more general than the RLPSI and is an optimum PPG-LPSI. In addition, it includes the LPSI and the RLPSI as particular cases.

Instead of using restriction $M'b=0$ to minimize $E[(H-1)^2]$, we can use restriction $C'b=\theta d$ and minimize $\Phi(b,v)=b'Pb+w'Gw-2w'Gb+2v'(C'b-\theta d)$ with respect to $b$, $v'$ and $\theta$ (Tallis 1985; Lin 2005) assuming that $P$, $G$, $U'$, $d$ and $w$ are known. The derivative results in matrix notation are

\begin{equation}\left( \begin{array}{c}
b_T\\
v \\
\theta 
\end{array} \right)=
\left( \begin{array}{ccc}
P & C & 0_{tx1} \\
C' & 0_{rxt} & -d \\
0'_{1xt} & -d' & 0
\end{array} \right)^{-1}
\left( \begin{array}{c}
Gw\\
\textbf{0} \\
0
\end{array} \right)\end{equation}

where $0_{tx1}$ is a null vector $tx1$, $0_{rxt}$ is a null matrix $rxt$ and $\textbf{0}$ is a null column vector $(r-1)x1$; $0$ is the standard zero value. The inverse matrix of coefficients 
$\left( \begin{array}{ccc}
P & C & 0_{tx1} \\
C' & 0_{rxt} & -d \\
0'_{1xt} & -d' & 0
\end{array} \right)^{-1}$ 

in Equation (4) is not easy to obtain; for this reason, Tallis (1985) obtained his results in two steps. That is, Tallis (1985) first derived Equation (3) with respect to $b$ and $v'$, from where he obtained $b_T=b_R+\theta \delta$ where $b_R=Kb$, $\delta=P^{-1}C(C'P^{-1}C)^{-1}d$ and $d=(d_1, d_2,..., d_r)$. Next he derived $E[(b'_Ty-H)^2]$ only with respect to $\theta$, and his result was $\theta=\frac{b'C(C'P^{-1}C)^{-1}d}{d'(C'P^{-1}C)^{-1}d}$ where $b=P^{-1}Gw$ is the LPSI vector of coefficients, C'=U'G, $d$ is the vector of the predetermined proportional gains imposed by the breeder and $P^{-1}$ is the inverse of matrix $P$. When $\theta=0$, $b_T=b_R$, and if $\theta=0$ and $U'$ is the null matrix, $b_T=b$. That is, the PPG-LPSI obtained by Tallis (1985) is more general that the RLPSI and the LPSI.

When $\theta=0$, Equation (5) is equal to $b_{T_0}=b_R+\delta$ the latter equation was the original result obtained by Tallis (1962). Tallis (1962) derived Equation (3) with respect to vectors $b$ and $v$ under the restriction $U'Gb=d$, i.e., without $\theta$ or $\theta=1$ . Later, James (1968) maximized the correlation between $I$ and $H$ ($\rho_{HI}$) under the Tallis (1962) restriction and once more obtained Equation (7). Mallard (1972) showed that Equation (7) is not optimum, i.e., it does not minimize $E[(I-H)^2]$ and does not maximize $\rho_{HI}$, and gave the optimum solution which we have presented here in Equation (2). Later, using restriction $U'Gb=\theta d$, Tallis (1985) obtained Equation (5), which also is optimum. 
	
Let $b_P=b_M=b_T$ be the PPG-LPSI vector of coefficients. Then, the optimum PPG-LPSI can be written as $I_P=b'_{P}y$ while the maximized correlation between the PPG-LPSI and the net genetic merit will be $\rho_{HI_P}=\frac{w'Gb_P}{\sqrt{w'Gw} \sqrt{b'_P P b_p}}$.

According to the conditions for constructing a valid PPG-LPSI, the index $I_P=b'_{P}y$ should have normal distributions. 

Under the predetermined restrictions imposed by the breeder, the $I_P=b'_{P}y$ should have maximum correlation with $H=w'g$ and it should be useful for ranking and selecting among individuals with different net genetic merit. However, for more than 2 restrictions the proportionality constant ($\theta$) could be lower than 1; in that case, $\rho_{HI_P}$ will be lower than the correlation between LPSI and $H=w'g$ ($\rho_{HI}$), in 
addition, when the restriction $M'b=0$ or $U'Gb=\theta d$ are imposed on the PPG-LPSI vector of coefficients, the restricted traits decrease their effect on the correlation between PPG-LPSI and $H=w'g$. 
	
The maximized PPG-LPSI selection response and expected genetic gains per trait can be written as $R_P=k_I \sqrt{b'_M P b_M}=k_I \sqrt{b'_T P b_T}$ and $E_P=k_I \frac{Gb_M}{\sqrt{b'_M P b_M}}=k_I \frac{Gb_T}{\sqrt{b'_T P b_T}}$ respectively, where $k_I$ is the standardized selection differential or selection intensity associated with the PPG-LPSI.

The maximized PPG-LPS selection response (Equation 10) has the same form as the maximized LPSI selection response. Thus, under   predetermined restrictions, Equation (10) predicts the mean improvement in $H$ due to indirect selection on $I_P=b'_P y$. Predetermined restriction effects will be observed on the PPG-LPSI expected genetic gain per trait (Equation 11).

### Statistical properties of the PPG-LPSI
Assuming that $H=w'g$ and $I_P=b'_P y$ have a bivariate joint normal distribution, $b_P=K_M b$, $b=P^{-1}Gw$, and $P$, $G$ and $w$ are known, the PPG-LPSI has the same properties as the RLPSI. Some of the main PPG-LPSI properties are:

\begin{itemize}
\item[1.]Matrices $Q_M=P^{-1}M(M'P^{-1}M)^{-1}M'$ and $K_M=[I-Q_M]$ have the same function as matrices $Q=P^{-1}C(C'P^{-1}C)^{-1}C'$ and $K=[I-Q]$ in the RLPSI.
\item[2.] Matrices $Q_M$ and $K_M$ are both projectors, i.e., they are idempotent ($K_M=K^2_M$ and $Q_M=Q^2_M$), unique and orthogonal, i.e., $K_M Q_M=Q_M K_M=0$.
\item[3.] Matrix $Q_M$ projects $b$ into a space generated by the columns of matrix $M$ due to the restriction $M'b=0$ that is introduced when $\Phi(b,v)$ is maximized with respect to $b$, while matrix $K_M$ projects $b$ into a space that is perpendicular to the space generated by the columns of matrix $M$ (Rao 2002). Thus the function of matrix $K_M$ is to transform vector $b=P^{-1}Gw$ into vector $b_P=K_M b$.
\item[4.] The variance of $I_P=b'_P y$ ($\sigma^2_{I_P}=b'_P P b_P$) is equal to the covariance between $I_P=b'_P y$ and $H=w'a$ ($\sigma_{HI_P}=w'G b_P$). As $K_M=K^2_M$, $K'_MP=PK_M$ and $b'P=w'G$, then, $\sigma^2_{I_P}=b'_P P b_P=b'K'_M b=b'PK^2_M b=w'Gb_P=\sigma_{HI_P}$.
\item[5.] The maximized correlation between $H$ and $I_P=b'_P y$ is equal to $\rho_{HI_P}=\frac{\sigma_{I_P}}{\sigma_H}$. In point 4 of this subsection, we showed that $\sigma_{HI_P}=\sigma^2_{I_P}$, then $\rho_{HI_P}=\frac{w'Gb_P}{\sqrt{w'Gw}\sqrt{b'_P P b_p}}=\sqrt{\frac{b'_P P b_p}{w'Gw}}=\frac{\sigma_{I_P}}{\sigma_H}$.
\item[6.] The variance of the predicted error, $Var(H-I_P)=(1-\rho^2_{HI_P})\sigma^2_H$, is minimal. By point 4 of this subsection, $\sigma_{HI_P}=\sigma^2_{I_P}$, then $Var(H-I_R)=\sigma^2_H-\sigma^2_{I_P}=(1-\rho^2_HI_P})\sigma^2_H$.
\item[7.] The heritability of the PPG-LPSI is equal to $h^2_{I_P}=\frac{b'_PG b_p}{b'_P P b_p}$.
\end{itemize}

### There only one optimum PPG-LGSI
Let $S=C'P^{-1}C$, under the restriction $D'd=0$, Itoh and Yamada (1987) showed that $D(D'SD)^{-1}D'=S^{-1}-S^{-1}d(d'S^{-1}d)^{-1}d'S^{-1}$, from where substituting $S^{-1}-S^{-1}d(d'S^{-1}d)^{-1}d'S^{-1}$ for $D(D'SD)^{-1}D'$ in matrix $Q_M$, Equation (2) can be written as Equation (5), i.e., $b_M=b_T$. Therefore, the Mallard (1972) and Tallis (1985) vectors of coefficients are the same. In addition, Itoh and Yamada (1987) showed that the Harville (1975) vector of coefficients can written as $\frac{b_T}{\sigma_{I_T}}$, where $\sigma_{I_T}$ is the standard deviation of the variance of Tallis (1985) PPG-LPSI. Thus, in reality there is only one optimum PPG-LPSI.

Itoh and Yamada (1987) also pointed out that matrix $D'=\left( \begin{array}{ccccc}
d_r & 0 & \cdots & 0 & -d_1 \\
0 & d_r & \cdots & 0 & -d_2 \\
\vdots & \vdots & \ddots & \vdots & \vdots\\
0 & 0 & \cdots & d_r & -d_{r-1} \end{array} \right)$
is only one example of several possible Mallard (1972) $D'$ matrices. They showed that any matrix $D'$ that satisfies condition $D'd=0$ is another Mallard (1972) matrix of predetermined proportional gains. According to Itoh and Yamada (1987), matrices:
\[D'=\left( \begin{array}{cccccc}
d_2 & -d_1 & 0 & \cdots & 0 & 0 \\
0 & d_3 & -d_2 & \cdots & 0 & 0 \\
\vdots & \vdots & \vdots & \ddots & \vdots & \vdots\\
0 & 0 & 0 & \cdots & d_r & -d_{r-1} \end{array} \right)\] 
\[D'=\left( \begin{array}{ccccc}
d_2 & -d_1 & 0 & \cdots & 0 \\
d_3 & 0 & -d_1 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots & \vdots\\
d_r & 0 & \cdots & 0 & -d_1 \end{array} \right)\] 
are also Mallard (1972) matrices of predetermined proportional gains because they satisfy condition $D'd=0$ . However, matrix $D'=\left( \begin{array}{ccccc}
d_r & 0 & \cdots & 0 & -d_1 \\
0 & d_r & \cdots & 0 & -d_2 \\
\vdots & \vdots & \ddots & \vdots & \vdots\\
0 & 0 & \cdots & d_r & -d_{r-1} \end{array} \right)$ is ''easier'' to construct.

Harville (1975) maximized the correlation between $I$ and $H$ ($\rho_{IH}$) under the restriction $C'b=\theta d$; he was the first to point out the importance of the proportionality constant ($\theta$) in the PPG-LPSI. Itoh and Yamada (1987) pointed out several problems associated with the Tallis (1985) PPG-PSI: 

\begin{itemize}
\item[(1)] When the number of restrictions imposed on the PPG-PSI expected genetic gains increases,$\theta$ tends to zero and then PPG-PSI accuracy decreases.
\item[(2)] The $\theta$ values could be negative, in which case PPG-PSI results have no meaning in practice.
\item[(3)] The PPG-PSI may cause the population means to shift in the opposite direction to the predetermined desired direction; this may happen due to the opposite directions between the economic values and the predetermined desired direction. Itoh and Yamada (1987) thought that one possible solution to those problems could be to use the linear phenotypic selection index with desired gains.
\end{itemize}

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

LPSI(datos,file.wgt,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",rawdata,one.env,block.ex,softR,file.covG)

```
[Return to examples](https://github.com/RAngelaPG/RIndSel-R/blob/master/Readme.md)
