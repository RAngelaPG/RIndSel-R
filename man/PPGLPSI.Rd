\name{PPGLPSI}
\Rdversion{1.0}
\alias{PPGLPSI}
\title{Predeterminated Proportional Gain Linear Phenotypic Selection Index.}
\description{
This function estimates the PPGLPSI.
}
\usage{
PPGLPSI(datos,file.wgt,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",
rawdata,one.env,block.ex,softR,file.covG)
}

\arguments{
  \item{datos}{dataframe object which contains the traits (columns) and entries (rows)}
  \item{file.wgt}{Text file in comma separated values (*.csv) format whose columns are: names of the traits, indicator variable for selected traits, economic weights or signs of the traits and indicator variable for restricted traits.}
  \item{selval}{Percent of selected genotypes (5, 10, 20, etc.). Default is selval=5.} 
  \item{design}{Field design ("lattice", "rcb", "AdjMeans"). Default is design="lattice".}
  \item{corr}{Selection index is based on variances-covariances matrix (FALSE) or correlation matrix (TRUE). Default is corr=FALSE.}
  \item{out}{Output file name. Default is out="out.txt".}
  \item{outcsv}{Output file in format *.csv. Default is outcsv="out.csv".}
  \item{rawdata}{If genotypic and phenotypic variances are already calculated, then it can be used to calculate directly the index by setting rawdata=FALSE. Default is TRUE.}
  \item{one.env}{Use FALSE for multienviromrent trials.}
  \item{block.ex}{Use FALSE always.}
  \item{softR}{Use "" always.}
  \item{file.covG}{When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.}
}

\details{
This index is called predetermined proportional gain linear phenotypic selection index (PPG-LPSI) because the breeder pre-sets optimal levels for certain traits before the selection is carried out. Some of the main objectives of the PPG-LPSI are
to optimize the expected genetic gain per trait, predict the net genetic merit \eqn{H=w'g} and select the individuals with the highest net genetic merit values as parents of the next generation. The PPG-LPSI allows imposing restrictions different 
from zero on the expected genetic gains of some traits, while other traits increase (or decrease) their expected genetic gains without imposing any restrictions. The PPG-LPSI solves the LPSI equations subject to the condition that the covariance 
between the LPSI and some linear functions of the genotypes involved be equal to a vector of predetermined constants or genetic gains defined by the breeder.

Let \eqn{d'=(d_1, d_2,...,d_r)} be a vector $rx1$ of the predetermined proportional gains and assume that \eqn{\mu_q} is the population mean of the \eqn{q^th} trait before selection. One objective could be to change \eqn{\mu_q} to \eqn{\mu_q+d_q}, 
where \eqn{d_q} is a predetermined change in \eqn{\mu_q} (in the RLPSI, \eqn{d_q = 0}, \eqn{q=1,2,...,r}, where \eqn{r} is the number of predetermined proportional gains). We can solve this problem in a similar manner as we did in the RLPSI. That is, 
minimizing the mean squared difference between \eqn{I} and \eqn{H} under the restriction \eqn{D'U'Gb = 0}, where \eqn{D'} 
\tabular{lllll}{
d_r \tab 0 \tab ... \tab 0 \tab -d_1 \cr
0 \tab d_r \tab ... \tab 0 \tab -d_2 \cr
. \tab . \tab . \tab . \tab . \cr
. \tab . \tab . \tab . \tab . \cr
. \tab . \tab . \tab . \tab . \cr
0 \tab 0 \tab ... \tab d_r \tab -d_{r-1}
}
is a Mallard (1972) matrix \eqn{(r-1)xr} of predetermined proportional gains, \eqn{d_q} (q=1,2,...,r) is the \eqn{q^th} element of vector \eqn{d'}, \eqn{U'} is the RLPSI matrix of restrictions of 1s and 0s, \eqn{G} is the covariance matrix of genotypic values and \eqn{b} is the LPSI vector of coefficients. Also, it is possible to minimize \eqn{E[(H-1)^2]} under the restriction \eqn{U'Gb = \theta d} (Tallis 1985), where \eqn{\theta} is a 
proportionality constant. Both approaches are very similar but the equations obtained when introducing the \eqn{D'U'Gb=0} restriction are simpler than when introducing \eqn{U'Gb = \theta d} restrictions into the process of minimizing \eqn{E[(H-1)^2]}.\cr

\bold{The maximized PPG-LPSI parameters}\cr 

Let \eqn{M'=D'C'} be the Mallard (1972) matrix of predetermined restrictions, where \eqn{C'=U'G}. Under the restriction \eqn{M'b=0}, we can minimize \eqn{E[(H-1)^2]} assuming that \eqn{P}, \eqn{G}, \eqn{U'}, \eqn{D'} and \eqn{w} are known; that is, we need to minimize the function \eqn{\Phi(b,v)=b'Pb+w'Gw-2w'Gb+2v'M'b} with respect to vectors \eqn{b} and \eqn{v'=(v_1, v_2,...,v_{r-1})}, where \eqn{v} is a vector of Lagrange multipliers. Equation (1) derivative results from \eqn{b} and \eqn{v}, i.e.,

\tabular{llllll}{
 P \tab M \tab | \tab b \tab \tab Gw\cr
  \tab  \tab | \tab  \tab = \tab \cr
 M' \tab 0 \tab | \tab v \tab \tab 0\cr
}

from where the vector that minimizes \eqn{E[(H-1)^2]} under the restriction \eqn{M'b = 0} is \eqn{b_M = K_Mb} where \eqn{K_M = [I_t-Q_M]}, \eqn{Q_M = P^-1 M (M'P^-1 M)^{-1} M' = P^-1 CD(D'C'P^-1 CD)^-1 D'C'} and \eqn{I_t} is an identity matrix of size \eqn{txt}. When \eqn{D = U}, \eqn{b_M = b_R} (the RLPSI vector of coefficients), and when \eqn{D=U} and \eqn{U'} is a null matrix, \eqn{b_M = b} (the LPSI vector of coefficients). Thus, the Mallard (1972) index is more general than the RLPSI and is an optimum PPG-LPSI. In addition, it includes the LPSI and the RLPSI as particular cases.

Instead of using restriction \eqn{M'b = 0} to minimize \eqn{E[(H-1)^2]}, we can use restriction \eqn{C'b = \theta d} and minimize
\eqn{\Phi(b,v) = b'Pb + w'Gw-2w'Gb + 2v'(C'b-\theta d)} with respect to \eqn{b}, \eqn{v'} and \eqn{\theta} (Tallis 1985; Lin 2005) assuming that \eqn{P}, \eqn{G}, \eqn{U'}, \eqn{d} and \eqn{w} are known. The derivative results in matrix notation are

\tabular{ccccccc}{
b_T\tab \tab P \tab C \tab 0_{tx1} \tab | \tab Gw \cr
v \tab =  \tab C' \tab 0_{rxt} \tab -d \tab | \tab 0\cr
\eqn{\theta} \tab \tab 0'_{1xt} \tab -d' \tab 0 \tab | \tab 0 \cr 
}

where \eqn{0_{tx1}} is a null vector \eqn{tx1}, \eqn{0_{rxt}} is a null matrix \eqn{rxt} and \eqn{0} is a null column vector \eqn{(r-1)x1}; \eqn{0} is the standard zero value. The inverse matrix of coefficients 

\tabular{ccc}{
P \tab C \tab 0_{tx1} \cr
C' \tab 0_{rxt} \tab -d \cr
0'_{1xt} \tab -d' \tab 0 \cr
}

in Equation (4) is not easy to obtain; for this reason, Tallis (1985) obtained his results in two steps. That is, Tallis (1985) first derived Equation (3) with respect to $b$ and \eqn{v'}, from where he obtained \eqn{b_T = b_R + \theta \delta} where \eqn{b_R = Kb}, \eqn{\delta = P^-1 C (C'P^-1 C)^-1 d} and \eqn{d=(d_1, d_2,..., d_r)}. Next he derived \eqn{E[(b'_Ty-H)^2]} only with respect to \eqn{\theta}, and his result was \eqn{\theta = (b'C(C'P^-1 C)^-1 d) / (d'(C'P^-1 C)^-1 d)} where \eqn{b = P^-1 Gw} is the LPSI vector of coefficients, \eqn{C'= U'G}, \eqn{d} is the vector of the predetermined proportional gains imposed by the breeder and \eqn{P^-1} is the inverse of matrix \eqn{P}. When \eqn{\theta=0}, \eqn{b_T=b_R}, and if \eqn{\theta = 0} and \eqn{U'} is the null matrix, \eqn{b_T = b}. That is, the PPG-LPSI obtained by Tallis (1985) is more general that the RLPSI and the LPSI.

When \eqn{\theta=0}, Equation (5) is equal to \eqn{b_T_0 = b_R + \delta} The latter equation was the original result obtained by Tallis (1962). Tallis (1962) derived Equation (3) with respect to vectors \eqn{b} and \eqn{v} under the restriction \eqn{U'Gb = d}, i.e., without \eqn{\theta} or \eqn{\theta = 1} . Later, James (1968) maximized the correlation between \eqn{I} and \eqn{H}(\eqn{\rho_HI}) under the Tallis (1962) restriction and once more obtained Equation (7). Mallard (1972) showed that Equation (7) is not optimum, i.e., it does not minimize \eqn{E[(I-H)^2]} and does not maximize \eqn{\rho_{HI}}, and gave the optimum solution which we have presented here in Equation (2). Later, using restriction \eqn{U'Gb=\theta d}, Tallis (1985) obtained Equation (5), which also is optimum. 
	
Let \eqn{b_P = b_M = b_T} be the PPG-LPSI vector of coefficients. Then, the optimum PPG-LPSI can be written as \eqn{I_P = b'_Py}
while the maximized correlation between the PPG-LPSI and the net genetic merit will be \eqn{\rho_HI_P= (w'Gb_P) \ (\sqrt w'Gw \sqrt b'_P P b_p)}

According to the conditions for constructing a valid PPG-LPSI, the index \eqn{I_P = b'_Py} should have normal distributions. 

Under the predetermined restrictions imposed by the breeder, the \eqn{I_P = b'_Py} should have maximum correlation with \eqn{H = w'g} and it should be useful for ranking and selecting among individuals with different net genetic merit. However, for more than 2 restrictions the proportionality constant (\eqn{\theta}) could be lower than 1; in that case, \eqn{\rho_HI_P} will be lower than the correlation between LPSI and \eqn{H=w'g} (\eqn{\rho_HI}), in addition, when the restriction \eqn{M'b = 0} or \eqn{U'Gb = \theta d} are imposed on the PPG-LPSI vector of coefficients, the restricted traits decrease their effect on the correlation between PPG-LPSI and \eqn{H = w'g}. 
	
The maximized PPG-LPSI selection response and expected genetic gains per trait can be written as
\deqn{R_P = k_I \sqrt b'_M P b_M = k_I \sqrt b'_T P b_T }
and
\deqn{E_P = k_I Gb_M / \sqrt b'_M P b_M = k_I Gb_T / \sqrt b'_T P b_T }
respectively, where \eqn{k_I} is the standardized selection differential or selection intensity associated with the PPG-LPSI.

The maximized PPG-LPS selection response (Equation 10) has the same form as the maximized LPSI selection response. Thus, under   predetermined restrictions, Equation (10) predicts the mean improvement in \eqn{H} due to indirect selection on 
\eqn{I_P = b'_P y}. Predetermined restriction effects will be observed on the PPG-LPSI expected genetic gain per trait (Equation 11).

\bold{Statistical properties of the PPG-LPSI}\cr

Assuming that \eqn{H = w'g} and \eqn{I_P = b'_P y} have a bivariate joint normal distribution, \eqn{b_P = K_M b}, \eqn{b = P^-1 Gw}, and \eqn{P}, \eqn{G} and \eqn{w} are known, the PPG-LPSI has the same properties as the RLPSI. Some of the main PPG-LPSI 
properties are:

\enumerate{
\item Matrices \eqn{Q_M = P^-1 M (M'P^-1 M)^-1 M'} and \eqn{K_M = [I-Q_M]} have the same function as matrices \eqn{Q = P^-1C(C'P^-1C)^-1 C'} and \eqn{K = [I-Q]} in the RLPSI.
\item Matrices \eqn{Q_M} and \eqn{K_M} are both projectors, i.e., they are idempotent (\eqn{K_M = K^2_M} and \eqn{Q_M = Q^2_M}), unique and orthogonal, i.e., \eqn{K_M Q_M = Q_M K_M = 0}.
\item Matrix \eqn{Q_M} projects \eqn{b} into a space generated by the columns of matrix \eqn{M} due to the restriction \eqn{M'b = 0} that is introduced when \eqn{\Phi(b,v)} is maximized with respect to \eqn{b}, while matrix \eqn{K_M} projects 
\eqn{b} into a space that is perpendicular to the space generated by the columns of matrix \eqn{M} (Rao 2002). Thus the function of matrix \eqn{K_M} is to transform vector \eqn{b = P^-1 Gw} into vector \eqn{b_P = K_M b}.
\item The variance of \eqn{I_P = b'_P y} (\eqn{\sigma^2_{I_P} = b'_P P b_P}) is equal to the covariance between \eqn{I_P = b'_P y} and \eqn{H=w'a} (\eqn{\sigma_{HI_P}=w'G b_P}). As \eqn{K_M=K^2_M}, \eqn{K'_MP = PK_M} and \eqn{b'P = w'G}, then, 
\eqn{\sigma^2_{I_P} = b'_P P b_P = b'K'_M b = b'PK^2_M b = w'Gb_P = \sigma_{HI_P}}.
\item The maximized correlation between \eqn{H} and \eqn{I_P = b'_P y} is equal to \eqn{\rho_HI_P = \sigma_{I_P} / \sigma_H}. In point 4 of this subsection, we showed that \eqn{\sigma_{HI_P} = \sigma^2_{I_P}}, 
then \eqn{\rho_{HI_P}= w'Gb_P / \sqrt w'Gw \sqrt b'_P P b_p = \sqrt{b'_P P b_p / w'Gw} = \sigma_{I_P} / \sigma_H}.
\item The variance of the predicted error, \eqn{Var(H-I_P) = (1-\rho^2_{HI_P})\sigma^2_H}, is minimal. By point 4 of this subsection, \eqn{\sigma_{HI_P} = \sigma^2_{I_P}}, then \eqn{Var(H-I_R) = \sigma^2_H-\sigma^2_{I_P} = (1-\rho^2_HI_P)\sigma^2_H}.
\item The heritability of the PPG-LPSI is equal to \eqn{h^2_{I_P} = b'_PG b_p / b'_P P b_p}.
}

\bold{There only one optimum PPG-LGSI}\cr

Let \eqn{S = C'P^{-1}C}, under the restriction \eqn{D'd = 0}, Itoh and Yamada (1987) showed that \eqn{D(D'SD)^{-1}D' = S^{-1}-S^{-1}d(d'S^{-1}d)^{-1}d'S^{-1}}, from where substituting \eqn{S^{-1}-S^{-1}d(d'S^{-1}d)^{-1}d'S^{-1}} for \eqn{D(D'SD)^{-1}D'} in matrix \eqn{Q_M}, Equation (2) can be written as Equation (5), i.e., \eqn{b_M=b_T}. Therefore, the Mallard (1972) and Tallis (1985) vectors of coefficients are the same. In addition, Itoh and Yamada (1987) showed that the Harville (1975) vector of coefficients can written as \eqn{b_T / \sigma_{I_T}}, where \eqn{\sigma_{I_T}} is the standard deviation of the variance of Tallis (1985) PPG-LPSI. Thus, in reality there is only one optimum PPG-LPSI.

Itoh and Yamada (1987) also pointed out that matrix \eqn{D'}:

\tabular{ccccc}{
d_r \tab 0 \tab ... \tab 0 \tab -d_1 \cr
0 \tab d_r \tab ... \tab 0 \tab -d_2 \cr
. \tab . \tab ... \tab . \tab . \cr
. \tab . \tab ... \tab . \tab . \cr
. \tab . \tab ... \tab . \tab . \cr
0 \tab 0 \tab ... \tab d_r \tab -d_{r-1} \cr
}

is only one example of several possible Mallard (1972) \eqn{D'} matrices. They showed that any matrix \eqn{D'} that satisfies condition \eqn{D'd=0} is another Mallard (1972) matrix of predetermined proportional gains. According to Itoh and 
Yamada (1987), matrices 
\eqn{D'=}
\tabular{cccccc}{
d_2 \tab -d_1 \tab 0 \tab ... \tab 0 \tab 0 \cr
0 \tab d_3 \tab -d_2 \tab ... \tab 0 \tab 0 \cr
. \tab . \tab . \tab ... \tab . \tab . \cr
. \tab . \tab . \tab ... \tab . \tab . \cr
. \tab . \tab . \tab ... \tab . \tab . \cr
0 \tab 0 \tab  0 \tab ... \tab d_r \tab -d_{r-1} \cr
}

\eqn{D'=}
\tabular{ccccc}{
d_2 \tab -d_1 \tab 0 \tab ... \tab 0 \cr
d_3 \tab 0 \tab -d_1 \tab ... \tab 0 \cr
. \tab . \tab ... \tab . \tab . \cr
. \tab . \tab ... \tab . \tab . \cr
. \tab . \tab ... \tab . \tab . \cr
d_r \tab 0 \tab ... \tab 0 \tab -d_1 \cr
}

are also Mallard (1972) matrices of predetermined proportional gains because they satisfy condition \eqn{D'd=0} . However, matrix \eqn{D'}

\tabular{ccccc}{
d_r \tab 0 \tab ... \tab 0 \tab -d_1 \cr
0 \tab d_r \tab ... \tab 0 \tab -d_2 \cr
. \tab . \tab ... \tab . \tab . \cr
. \tab . \tab ... \tab . \tab . \cr
. \tab . \tab ... \tab . \tab . \cr
0 \tab 0 \tab ... \tab d_r \tab -d_{r-1} \cr
}
is ''easier'' to construct.

Harville (1975) maximized the correlation between \eqn{I} and \eqn{H} (\eqn{\rho_IH}) under the restriction \eqn{C'b = \theta d}; he was the first to point out the importance of the proportionality constant (\eqn{\theta}) in the PPG-LPSI. Itoh and Yamada (1987) pointed out several problems associated with the Tallis (1985) PPG-PSI: 

\enumerate{
\item When the number of restrictions imposed on the PPG-PSI expected genetic gains increases,\eqn{\theta} tends to zero and then PPG-PSI accuracy decreases.
\item The \eqn{\theta} values could be negative, in which case PPG-PSI results have no meaning in practice.
\item The PPG-PSI may cause the population means to shift in the opposite direction to the predetermined desired direction; this may happen due to the opposite directions between the economic values and the predetermined desired direction. 
Itoh and Yamada (1987) thought that one possible solution to those problems could be to use the linear phenotypic selection index with desired gains.
}

}

\value{
The PPG-LPSI selection index and an output file which contains the selected entries and other useful information.
}

\examples{
library(Rindsel)
datos<-data.frame(read.csv("https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/C1_PSI_05_Phen.csv",header=T,na.strings=c(NA,".","-"))) #Raw data to analized.
datos$REP=as.factor(datos$REP)                                                                     #transform to factor.
datos$Block=as.factor(datos$Block)                                                                 #transform to factor.
datos$ENTRY=as.factor(datos$ENTRY)                                                                 #transform to factor.
file.wgt<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/weights_C1_PSI.csv"     #name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                          #Selection intensity.
design<-"lattice"                                                                                  #Experimental design.
corr<-FALSE                                                                                        #You can decide if you want to work ith the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                      #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
one.env<-TRUE                                                                                      #Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                                    #Use FALSE always.
softR<-""                                                                                          #Use "" always.
file.covG<-""                                                                                      #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.
PPGLPSI(datos,file.wgt,selval,design,corr,out="outextPPGLPSI.txt",outcsv="outPPGLPSI.csv",rawdata,one.env,block.ex,softR,file.covG)
file.show("outextPPGLPSI.txt")
}


\author{
Maintainer: Angela Pacheco <r.a.pacheco@cgiar.org>
Sergio Perez-Elizalde <sergiop@colpos.mx>
Gregorio Alvarado,
Jesus Ceron,
Francisco Rodriguez,
Jose Crossa <j.crossa@cgiar.org>
Juan Burgueno <j.burgueno@cgiar.org>
}

\references{
\enumerate{
\item Ceron-Rojas, J.J., J. SahagÃºn-Castellanos, F. Castillo-Gonzalez, and A. Santacruz-Varela, 2006. A selection index method based on eigenanalysis. Crop Science 46: 1711-1721.\cr
\item Ceron-Rojas, J.J., J. SahagÃºn-Castellanos, F. Castillo-Gonzalez, A. Santacruz-Varela and J. Crossa, 2008. A restricted selection index method based on eigenanalysis. Journal of Agricultural, Biological, and Environmental Statistics 13: 449-457.\cr
\item Cunningham, M. D., Mace, J. W., and Peters, E. R. (1969).Clinical experience with phenobarbitone in icterus neonatorum. Lancet, 1, 550.\cr
\item Bulmer, M. G. 1986. Sex ratio theory in geographically structured populations. Heredity, 56, 69-73.\cr
\item Hazel, L.N. 1943 The genetic basis for constructing selection indexes. Genetics 28, 476.\cr
\item Jansen RC. Studying complex biological systems using multifactorial perturbation. Nat Rev Genet 2003, 4: 145-151.\cr
\item Kempthorne, O., and A. W. Nordskog, 1959. Restricted selection indices. Biometrics 15: 10-19.\cr
\item Lande, R., and R. Thompson, 1990. Efficiency of marker-assisted in the improvement of quantitative traits. Genetics 124: 743-756.\cr
\item Li Z 1998. Molecular analysis of epistasis affecting complex traits. In: Paterson AH (ed) Molecular dissection of complex traits. CRC Press, Boca Raton, pp 119-130\cr
\item Muirhead, R. J., 2005. Aspects of Multivariate Statistical Theory, John Wiley and Sons, New York.\cr
\item Smith, H. F., 1936. A discriminant function for plant selection, In: Papers on Quantitative Genetics and Related Topics. Deparment of Genetics, North Carolina State College,Raleigh, North Carolina, pp. 466-476.
}
}

\keyword{ package }
