\name{BLPSI}
\Rdversion{1.0}
\alias{BLPSI}
\title{Base Linear Phenotypic Selection Index.}
\description{
This function estimates the BLPSI selection index.
}
\usage{
BLPSI(datos,file.wgt,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",
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
To derive the LPSI theory, we assumed that the phenotypic (\eqn{P}) and the genotypic (\eqn{G}) covariance matrix, and the vector of economic values (\eqn{w}) are known. However, \eqn{P},\eqn{G} and \eqn{w} are generally unknown and it is necessary to estimate them. There are many methods for estimating \eqn{P} and \eqn{G} (Lynch and Walsh 1998) and \eqn{w} (Cotterill and Jackson 1985; Magnussen 1990). However, when the estimator of \eqn{P} is not positive definite (all eigenvalues positive) or the estimator of G is not positive semidefinite (no negative eigenvalues), the estimator of G could be biased. In this case, the base linear phenotypic selection index (BLPSI): 
\deqn{I_B=w'y}
may be a better predictor of \eqn{H=w'g} than the estimated LPSI \eqn{I = b'y} (Williams 1962a; Lin 1978) if indeed the vector of economic values is known. Many authors (Williams 1962b; Harris 1964; Hayes and Hill 1980, 1981) have investigated the influence of parameter estimation errors on LPSI accuracy and concluded that those errors affect the accuracy of \eqn{I = b'y} when the accuracy of \eqn{P} and \eqn{G} is low. If vector \eqn{w} values are known, the BLPSI has certain advantages because of its simplicity and its freedom from parameter estimation errors (Lin 1978). Williams (1962a) pointed out that the BLPSI is superior to \eqn{I = \hatb'y} unless a large amount of data is available for estimating \eqn{P} and \eqn{G}.

There are some problems associated with the BLPSI. For example, what is the BLPSI selection response and the BLPSI expected genetic gains per trait when no data are available for estimating \eqn{P} and \eqn{G}? The BLPSI is a better selection index than the standard LPSI only if the correlation between the BLPSI and the net genetic merit is higher than the correlation between the LPSI and the net genetic merit (Hazel 1943). But if estimations of \eqn{P} and \eqn{G} are not available, how can we obtain the correlation between the base index and the net genetic merit? Williams (1962b) pointed out that the correlation between the BLPSI and \eqn{H=w'g} can be written as
\deqn{\rho_ HI_B = \sqrt w'Gw/w'Pw}
and indicated that the ratio \eqn{\rho_HI_B/\rho_HI} can be used to compare LPSI efficiency vs BLPSI efficiency; however, in the latter case we at least need to know the estimates of \eqn{P} and \eqn{G}, i.e., \eqn{\hat{P}} and \eqn{\hat{G}}.

In addition, Equation (1) is only an assumption, not a result, and implies that \eqn{P} and \eqn{G} are the same. That is, \eqn{b=P^{-1}Gw=w} only when \eqn{P=G}, which indicates that the BLPSI is a special case of the LPSI. Thus to obtain the selection response and the expected genetic gains per trait of the BLPSI, we need some information about \eqn{P} and \eqn{G}. Assuming that BLPSI is indeed a particular case of the LPSI, the BLPSI selection response and the BLPSI expected genetic gains per trait could be written as 
\deqn{R_B = k_I \sqrt w'Pw}
and 
\deqn{E_B = k_IGw / \sqrt w'Pw}
respectively. The parameters of Equations (3) and (4) were defined earlier.

There are additional implications if \eqn{b=P^-1Gw}. For example, if \eqn{P=G}, then \eqn{\rho_HI_B = \sqrt w'Gw / w'Pw} and BLPSI heritability \eqn{h^2_I_B = w'Gw / w'Pw} are equal to 1. However, in practice, the estimated values of the \eqn{\rho_HI_B} are usually lower than the estimated values of the \eqn{\rho_HI}.\cr
}

\value{
The BLPSI selection index and an output file which contains the selected entries and other useful information.
}

\examples{
library(Rindsel)
datos<-data.frame(read.csv("https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/C1_PSI_05_Phen.csv",header=T,na.strings=c(NA,".","-"))) #Raw data to analized.
datos$REP=as.factor(datos$REP)																		#transform to factor.
datos$Block=as.factor(datos$Block)																	#transform to factor.
datos$ENTRY=as.factor(datos$ENTRY)																	#transform to factor.
file.wgt<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/weights_C1_PSI.csv"     #name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                    		#Selection intensity.
design<-"lattice"                                                                            		#Experimental design.
corr<-FALSE                                                                                  		#You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                		#By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
one.env<-TRUE                                                                                		#Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                              		#Use FALSE always.
softR<-""                                                                                    		#Use "" always.
file.covG<-""                                                                                		#When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.
BLPSI(datos,file.wgt,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",rawdata,one.env,block.ex,softR,file.covG)
file.show("outextBLPSI.txt")
}


\author{
Maintainer: Angela Pacheco <r.a.pacheco@cgiar.org>
Sergio Perez-Elizalde <sergiop@colpos.mx>
Gregorio Alvarado
Jesus Ceron
Francisco Rodriguez
Jose Crossa <j.crossa@cgiar.org>
Juan Burgueno <j.burgueno@cgiar.org>
}

\references{
\enumerate{
\item Ceron-Rojas, J.J., J. Sahagún-Castellanos, F. Castillo-Gonzalez, and A. Santacruz-Varela, 2006. A selection index method based on eigenanalysis. Crop Science 46: 1711-1721.\cr
\item Ceron-Rojas, J.J., J. Sahagún-Castellanos, F. Castillo-Gonzalez, A. Santacruz-Varela and J. Crossa, 2008. A restricted selection index method based on eigenanalysis. Journal of Agricultural, Biological, and Environmental Statistics 13: 449-457.\cr
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
