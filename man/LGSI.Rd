\name{LGSI}
\Rdversion{1.0}
\alias{LGSI}
\title{Linear Genomic Selection Index.}
\description{
This function estimates the LGSI selection index.
}
\usage{
LGSI(datos,file.wgt,selval,design,corr,method,out="outextLGSI.txt",outcsv="outLGSI.csv",
rawdata,file_nameMARK,one.env,block.ex,softR,file.covG,file_nameTMARK,LG)
}

\arguments{
  \item{datos}{dataframe object which contains the traits (columns) and entries (rows)}
  \item{file.wgt}{Text file in comma separated values (*.csv) format whose columns are: names of the traits, indicator variable for selected traits, economic weights or signs of the traits and indicator variable for restricted traits.}
  \item{selval}{Percent of selected genotypes (5, 10, 20, etc.). Default is selval=5.} 
  \item{design}{Field design ("lattice", "rcb", "AdjMeans"). Default is design="lattice".}
  \item{corr}{Selection index is based on variances-covariances matrix (FALSE) or correlation matrix (TRUE). Default is corr=FALSE.}
  \item{method}{Choose between "vanraden" or "classic" for create the relationship marker matrix.}
  \item{out}{Output file name. Default is out="out.txt".}
  \item{outcsv}{Output file in format *.csv. Default is outcsv="out.csv".}
  \item{rawdata}{If genotypic and phenotypic variances are already calculated, then it can be used to calculate directly the index by setting rawdata=FALSE. Default is TRUE.}
  \item{file_nameMARK}{Text file in comma separated values (*.csv) format where we can find training marker information.}
  \item{one.env}{Use FALSE for multienvironment trials.}
  \item{block.ex}{Use FALSE always.}
  \item{softR}{Use "" always.}
  \item{file.covG}{When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.}
  \item{file_nameTMARK}{Text file in comma separated values (*.csv) format where we can find testing marker information.}
  \item{LG}{Interval between selection cycles}
}

\details{
In a similar manner as the linear phenotypic selection index (LPSI), the objective of the LGSI is to predict the net genetic merit \eqn{H=w'g}, where \eqn{g'=(g_1,g_2,...,g_t)} (t number of traits) is a vector of unobservable 
true breeding values and \eqn{w'=(w_1,w_2,...,w_t)} is a vector of economic weights. Suppose that the genomic breeding values \eqn{\gamma_i=Xu_i} are known; then, the LGSI can be written as
\eqn{I_G=\beta' \gamma}
where \eqn{\beta} is an unknown vector of weights.

The main advantage of the LGSI over the LPSI lies in the possibility of reducing the intervals between selection cycles (\eqn{L_G}) by more than two thirds (Lorenz et al. 2011); thus, this parameter should be incorporated into the LGSI selection
response and the expected genetic gain per trait to reflect the main advantage of the LGSI over the LPSI and the other indices. In the LPSI context we wrote the selection response as \eqn{R_I=k_I\sigma_H \rho_{HI}}, while the LGSI selection 
response can be written as
\eqn{R_{I_G}=\frac{k_I \sigma_{HI_G}}{L_G \sigma^2_{I_G}}=\frac{k_I}{L_G}\sigma_H \rho_{HI_G}}
where \eqn{k_I} is the standardized selection differential (or selection intensity) associated with LGSI, \eqn{\sigma_{HI_G}} is the covariance between \eqn{H=w'g} and LGSI, \eqn{\sigma^2_{I_G}} is the variance of LGSI, \eqn{\sigma_H} is the 
standard deviation of \eqn{H}, \eqn{\rho_{HI_G}} is the correlation between \eqn{H} and LGSI, and \eqn{L_G} denotes the intervals between selection cycles. 
	
Let \eqn{C} and \eqn{\Gamma} be matrices of covariance of the breeding values (\eqn{g}) and of the genomic breeding values (\eqn{\gamma}), respectively; then, the correlation between \eqn{H=w'g} and \eqn{I_G=\beta' \gamma} can be written as
\eqn{\rho_{HI_G}=\frac{w' \Gamma \beta}{\sqrt{w'Cw}\sqrt{\beta' \Gamma \beta}}}
where \eqn{w' \Gamma \beta= \sigma_{HI_G}} is the covariance between \eqn{H=w'g} and \eqn{I_G=\beta' \gamma} , \eqn{\sigma_H=\sqrt{w'Cw}} is the standard deviation of the variance of \eqn{H=w'g}, and \eqn{\sigma_{I_G}=\sqrt{\beta' \Gamma \beta}}
is the standard deviation of the variance of \eqn{I_G=\beta' \gamma}.

\bold{The maximize LGSI parameters.}\cr

To maximize the genomic selection response (Equation 2), suppose that \eqn{k_I}, \eqn{\sigma_H} and \eqn{L_G} are fixed and take the derivative of the natural logarithm (ln) of the correlation between \eqn{H} and \eqn{I_G} (Equation 3) with 
respect to vector \eqn{\beta}, equate the result of the derivative to the null vector, and isolate \eqn{\beta}, i.e.,
\eqn{\frac{d}{d\beta}ln \rho_{HI_G}=\frac{d}{d\beta}ln \left(\frac{w'\Gamma \beta}{\sqrt{w'Cw}\sqrt{\beta' \Gamma \beta}} \right)=0}
The result is \eqn{\beta=sw}, where \eqn{s=\frac{\beta' \Gamma \beta}{w' \Gamma \beta}} is a proportional constant that does not affect the maximum value of \eqn{\rho_{HI_G}} because this is invariant to the scale change; then assuming that 
\eqn{\beta=w}, the maximized LGSI selection response can be written as
\eqn{R_{I_G}=\frac{k_I}{L_G}\sqrt{w' \Gamma w}}

Hereafter we will refer to the LGSI genomic selection response as that of Equation (5). 

Also, because \eqn{\beta=w}, Equation (3) can be written as
\eqn{\rho_{HI_G} = \sqrt (w' \Gamma w) / \sqrt (w' C w) = \sigma_{I_G} / \sigma_H}
which is the maximized correlation between \eqn{H = w'g} and \eqn{I_G = \beta' \gamma}, or LGSI accuracy; \eqn{\sigma_H= \sqrt w'Cw } is the standard deviation of the variance of \eqn{H}, and \eqn{\sigma_{I_G} = \sqrt \beta' \Gamma \beta } is the 
standard deviation of the variance of \eqn{I_G}.

The LGSI expected genetic gain per trait (\eqn{E_{I_G}}) can be written as
\eqn{E_{I_G}=\frac{k_I \Gamma w}{L_G \sqrt{w' \Gamma w}}}
All the terms in Equation (7) were previously defined.

Let \eqn{\lambda_G=\frac{\rho_{HI_G}}{\rho_{HI}}} be LGSI efficiency vs. LPSI efficiency to predict the net genetic merit, where \eqn{\rho_{HI_G}} is the LGSI accuracy and \eqn{\rho_{HI}} the LPSI accuracy; in percentage terms, LGSI efficiency 
vs LPSI efficiency for each selection cycle can be written as
\eqn{p_G=100(\lambda_G-1)}

According to Equation (8), if \eqn{p_G>0}, LGSI efficiency will be greater than LPSI efficiency; if \eqn{p_G=0}, the efficiency of both selection indices will be equal, and if \eqn{p_G<0}, LPSI will be more efficient than LGSI for predicting 
\eqn{H=w'g}.

Equation (8) is useful for measuring LGSI efficiency in terms of accuracy when predicting the net genetic merit (\eqn{H=w'g}), while the Technow et al. (2013) inequality measures LGSI efficiency in terms of the time needed to complete one 
selection cycle. In the context of LGSI and LPSI, the Technow inequality can be written as
\eqn{L_G<\frac{\rho_{HI_G}}{h_I}L_P}
where \eqn{L_G} and \eqn{L_P} denote the time required to complete one selection cycle for LGSI and LPSI, respectively; \eqn{\rho_{HI_G}} is the LGSI accuracy, and \eqn{h_I} is the square root of the heritability (Lin and Allaire 1977; 
Nordskog 1978) of LPSI, which can be denoted as \eqn{h_I=\sqrt{\frac{b'Cb}{b'Pb}}}. Then, assuming that the selection intensity is the same for both selection indices, if Equation (9) is true, LGSI will be more efficient than LPSI per unit 
of time.

\bold{Statistical LGSI properties.}\cr

Assuming that \eqn{H} and \eqn{I_G} have joint bivariate normal distribution and that \eqn{\Gamma}, \eqn{C} and \eqn{w} are known, the LGSI has the following properties:

\enumerate{
\item The variance of \eqn{I_G} (\eqn{\sigma^2_{I_G}}) and the covariance between \eqn{H} and \eqn{I_G} (\eqn{\sigma_{HI_G}}) are equal, i.e., \eqn{\sigma^2_{I_G}=\sigma_{HI_G}}.
\item The maximized correlation between \eqn{H} and \eqn{I_G} (or LGSI accuracy) is equal to \eqn{\rho_{HI_G}=\frac{\sigma_{I_G}}{\sigma_H}}, where \eqn{\sigma_{I_G}} is the standard deviation of \eqn{\sigma^2_{I_G}} and \eqn{\sigma_H}
is the standard deviation of the variance of \eqn{H} (\eqn{\sigma^2_H}).
\item The variance of the predicted error, \eqn{Var(H-I_G)=(1-\rho^2_{HI_G})\sigma^2_H}, is minimal. Note that \eqn{Var(H-I_G)=\sigma^2_{I_G}+\sigma^2_H-2\sigma_{HI_G}}, and when \eqn{\beta=w}, \eqn{\sigma^2_{I_G}=\sigma_{HI_G}}, from where 
\eqn{Var(H-I_G)=\sigma^2_H-\sigma^2_{I_G}=(1-\rho^2_{HI_G})\sigma^2_H} is minimal.
\item The total variance of \eqn{H} explained by \eqn{I_G} is \eqn{\sigma^2_{I_G}=\rho^2_{HI_G}\sigma^2_H}. It is evident that if \eqn{\rho_{HI_G}=1}, \eqn{\sigma^2_{I_G}=\sigma^2_H}, and if \eqn{\rho_{HI_G}=0}, \eqn{\sigma^2_{I_G}=0}. 
That is, the variance of \eqn{H} explained by \eqn{I_G} is proportional to \eqn{\rho_{HI_G}}, and when \eqn{\rho_{HI_G}} is close to 1, \eqn{\sigma^2_{I_G}} is close to \eqn{\sigma^2_{H}}; if \eqn{\rho_{HI_G}} is close to 0, \eqn{\sigma^2_{I_G}}
is close to 0.
}

The LGSI properties described in points 1 to 4 of this subsection are the same as the LPSI properties described. This corroborates that the LGSI is an application of the LPSI theory to the GS context.

}

\value{
The LGSI selection index and an output file which contains the selected entries and other useful information.
}



\examples{
library(Rindsel)
datos<-data.frame(read.csv("https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/Data_Phenotypes_LGSI.csv",header=T,na.strings=c(NA,".","-"))) #Raw data to analized.
file.wgt<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/weights_LGSI.csv"   			#name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                    			#Selection intensity.
design<-"lattice"                                                                            			#Experimental design.
corr<-FALSE                                                                                  			#You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
method<-"vanraden"																				        #For create the relationship marker matrix.
rawdata<-TRUE                                                                                			#By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
file_nameMARK<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/Training_population_LGSI.csv" #name of the file training markers information.
one.env<-TRUE                                                                                			#Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                              			#Use FALSE always.
softR<-""                                                                                    			#Use "" always.
file.covG<-""                                                                                			#When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.
file_nameTMARK<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/Testing_population_LGSI.csv"	#name of the file testing markers information.
LG<-1																									#Interval between selection cycles.
LGSI(datos,file.wgt,selval,design,corr,method,out="outextLGSI.txt",outcsv="outLGSI.csv",rawdata,file_nameMARK,one.env,block.ex,softR,file.covG,file_nameTMARK,LG)
file.show("outextLGSI.txt")
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
