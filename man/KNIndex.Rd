\name{KNIndex}
\Rdversion{1.0}
\alias{KNIndex}
\title{Restrictive Linear Phenotypic Selection Index}
\description{
Funtion for the restrictive Kempthorne & Nordskog selection index}
\usage{
KNIndex(datos,file.wgt,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",
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
Kempthorne and Nordskog (1959) proposed an index that incorporates 
restrictions based on a predeteined level and improvement 
developed the theory of restricted selection index 
(ISR), which maximizes the genetic progress of some average 
Phenotypic variables while the other half is maintained without 
exchange. According to Cunningham {et al. } (1970), it is considered that an IS is 
restricted if the average allows one or more variables 
phenotype remains unchanged (no genetic gain) while 
in other genetic gain is maximized average. In general, the 
Kempthorne and Nordskog ISR (1959) has similar properties to IS 
Smith (1936). 

Consider the problem of increasing grain yield in 
corn while other phenotypic variables remain without 
exchange. In this situation Kempthorne and Nordskog (1959) (KN) proposed 
maximize \eqn{\rho_{YZ}^2}  incorporating constraints matrix 
genotypic variance-covariance. Bulmer (1980) showed that 
Kempthorne and Nordskog results (1959) also obtained 
maximize \eqn{ {\theta}'G \beta } . In the derivation of the ISR 
Kempthorne and Nordskog (1959) follows the approach of Bulmer (1980). 

Suppose we have \eqn{q} variables and desired phenotypic 
improve the average of \eqn{qr} of them while \eqn{ r } of them remain without 
average exchange, the index \eqn{Y=\beta'p } must maximize 
 \eqn{\rho_{YZ} ^ 2 } under such conditions. The restrictions are introduced by 
through the matrix \eqn{ W } of 0s and 1s, where 1s indicate that 
mean values of the variables remain unchanged, while the 
0s indicate average changes in phenotypic variables, such 
so that \eqn{ W'G \beta =0 } . Let \eqn{  C =GW } , since \eqn{ \rho_{YZ} } is invariant under changes of scale, can 
assume that \eqn{ \beta 'P \beta = 1 } , so \eqn{ \theta 'G \beta } should be maximized under the 
constraint \eqn{\beta'C = 0} and \eqn{\beta' P \beta = 1 } . Let \eqn{ \tau } and 
 \eqn{ v'=\left [{v_1 ... v_r} \right] } Lagrange multipliers. 
Should be maximized 
\deqn{ 
\Psi =\theta'G \beta - \tau  (\beta'P \beta-1) -v 'C' \beta}  
In deriving \eqn{ \Psi } over \eqn{ \beta } and equate the result to 
null vector is 
\deqn{ 
\theta'G-2 \tau P \beta -v'C' =0 
} 
Thus, the vector of coefficients of KN that maximizes \eqn{ \rho_{YZ}^2 } is 
 \eqn{ \beta_{KN}=[I -P^{-1}C (C'P^{-1}C)^{-1}C']\beta_S } , where \eqn{P^{-1} } 
 is the inverse matrix \eqn{P} and \eqn{\beta_S } is the vector defined in Equation 3 (the 
subscript denotes Kempthorne and Nordskog KN). If \eqn{ A=[ I - P^{-1} C(C'P ^{-1}C)^{-1}C'] } , then 
\deqn{ 
 \beta_{KN} =A \beta_S.  
}
}
\value{
The KN selection index and an output file which contains the selected entries and other useful infoation.
}

\examples{
library(Rindsel)
datos<-data.frame(read.csv("https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/C1_PSI_05_Phen.csv",header=T,na.strings=c(NA,".","-"))) #Raw data to analized.
datos$REP=as.factor(datos$REP)                                                                     #transform to factor.
datos$Block=as.factor(datos$Block)                                                                 #transform to factor.
datos$ENTRY=as.factor(datos$ENTRY)                                                                 #transform to factor.
file.wgt<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/weights_C1_PSI.csv"     #name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                    		   #Selection intensity.
design<-"lattice"                                                                            		   #Experimental design.
corr<-FALSE                                                                                  		   #You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                		   #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
one.env<-TRUE                                                                                		   #Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                              		   #Use FALSE always.
softR<-""                                                                                    		   #Use "" always.
file.covG<-""                                                                                		   #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.
KNIndex(datos,file.wgt,selval,design,corr,out="outextKN.txt",outcsv="outKN.csv",rawdata,one.env,block.ex,softR,file.covG)
file.show("outextKN.txt")
}


\author{
Maintainer: Angela Pacheco <r.a.pacheco@cgiar.org>
Sergio PÈrez-Elizalde <sergiop@colpos.mx>
Gregorio Alvarado,
Jesus Ceron,
Francisco Rodriguez,
JosÈ Crossa <j.crossa@cgiar.org>
Juan Burgueno <j.burgueno@cgiar.org>
}

\references{
\enumerate{
\item Ceron-Rojas, J.J., J. Sahag√∫n-Castellanos, F. Castillo-Gonzalez, and A. Santacruz-Varela, 2006. A selection index method based on eigenanalysis. Crop Science 46: 1711-1721.\cr
\item Ceron-Rojas, J.J., J. Sahag√∫n-Castellanos, F. Castillo-Gonzalez, A. Santacruz-Varela and J. Crossa, 2008. A restricted selection index method based on eigenanalysis. Journal of Agricultural, Biological, and Environmental Statistics 13: 449-457.\cr
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
