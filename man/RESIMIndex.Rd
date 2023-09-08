\name{RESIMIndex}
\Rdversion{1.0}
\alias{RESIMIndex}
\title{Restrictive Eigen Selection Index Method}
\description{
Funtion for restrictive eigen selection index method
}
\usage{
RESIMIndex(datos,file.wgt,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",
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
  \item{file.covG}{When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.}}

\details{
RESIM is a generalization of ESIM to the case in which are incorporated 
restrictions on \eqn{Y ={\rm{\bf{\beta} 'p}} }{Y=\beta'p} , the procedure for 
set the \eqn{ {\rm{\bf \beta}} }{\beta} is, therefore, very similar to ESIM 
and is based on IS restricted Kempthorne and Nordskog (1959). 
So, suppose you have \eqn{q} variables and phenotypic 
want to improve the means of \eqn{q-r} of them while \eqn{r} of them 
average remain unchanged, the index \eqn{Y ={\rm{\bf{\beta}'p}} }{Y=\beta'p} 
should maximize \eqn{\rho_{YZ}^2} under conditions tale. Again, the 
restrictions are introduced by the matrix \eqn{ {\rm{\bf W}} }{W} of 0s and 
1s, where 1s indicate that the mean values of the variables 
Phenotypic remain unchanged, while indicating 0s 
average changes in phenotypic variables, such that \eqn{ {\rm 
{\bf{W}'G \beta}} ={\rm{\bf 0}} }{W'G \beta=0} . Let \eqn{ {\rm{\bf C}} ={\rm{\bf}} }{C=} GW, 
since \eqn{ \rho_{YZ} } is invariant under changes of scale can be assumed 
which \eqn{ {\rm{\bf{\beta}'P \beta}} = 1 }{\beta'*P*\beta=1} and \eqn{ {\rm{\bf{\theta}'G \theta 
}}=1 }{\theta'G\theta=1} , so \eqn{ \rho_{YZ}^2 } be maximized under the constraints 
 \eqn{ {\rm{\bf{\beta}'C}} ={\rm{\bf 0}} }{\beta'C=0} , \eqn{ {\rm{\bf{\beta}'P \beta}} = 1 }{\beta'*P*\beta=1} and 
 \eqn{ {\rm{\bf{\theta}'G \theta}}=1 }{\theta'G\theta=1} . Let \eqn{ \mu } , \eqn{ \omega } and \eqn{ {\rm{\bf 
{v}'}}=\left [{v_1 ... v_r} \right] }{v'=[v_1 ... v_r]} Lagrange multipliers. It should 
maximized 
\deqn{ 
\Phi = \left ({{\rm{\bf{\theta} 'G \beta}}} \right) ^ 2 - \mu \left ({{\rm{\bf 
{\beta} 'P \beta}} -1} \right) - \omega \left ({{\rm{\bf{\theta}' G \theta 
-1}}} \right) -{\rm{\bf{v} '{C}' \beta}} 
}{\Phi = (\theta '_M*K*\beta_M-1 )^2-\mu(\beta '_M*T*\beta_M-1) -\omega(\theta '_M*K*\theta_M-1 )-v'C'\beta} 
In deriving \eqn{ \Phi } over \eqn{ {\rm{\bf \beta}} }{\beta} , \eqn{ {\rm{\bf \theta 
}} }{\theta} , \eqn{ \mu } , \eqn{ \omega } , and \eqn{ {\rm{\bf{v} '}}=\left [{v_1 ... v_r} \right] }{v'=[v_1 ... v_r]} , 
Ceron-Rojas {et al}. (2008a) showed that 
\deqn{ 
({\rm{\bf Q}} - \omega{\rm{\bf I}}){\rm{\bf \beta}} ={\rm{\bf 0}} 
}{(Q-\omega I)\beta=0}  
where now \eqn{ {\rm{\bf Q}} = [{\rm{\bf I}} -{\rm{\bf P}}^{-1}{\rm{\bf 
C}} ({\rm{\bf{C}'P}}^{-1}{\rm{\bf C}}) ^{-1}{\rm{\bf{C}'}}]{\rm{\bf 
P}}^{-1}{\rm{\bf G}} }{Q=[I-P^{-1}C(C'P^{-1}C)^{-1}P^{-1}G]} . Thus, RESIM, the value that maximizes \eqn{ \rho
_{YZ}^2 } under the constraint \eqn{ {\rm{\bf{\beta}'P \beta}} = 1 }{\beta'P \beta = 1} , \eqn{ {\rm 
{\bf{\theta}'G \theta}} = 1 }{\theta'G \theta = 1} and \eqn{ {\rm{\bf{\beta}'C}} ={\rm{\bf 0}} }{\beta'C=0} is the 
first eigenvalue ( \eqn{ \omega) } of the matrix \eqn{ {\rm{\bf Q}} }{Q} , and the vector 
lets build resim \eqn{Y={\rm{\bf{\beta}'p}} }{Y=\beta'p} (with maximum 
correlation with \eqn{Z={\rm{\bf{\theta}'g}}) }{Z=\theta'g} is the first eigenvector 
( \eqn{ {\rm{\bf \beta}}) }{\beta} matrix \eqn{ {\rm{\bf Q}} }{Q}
}

\value{
The RESIM selection index and an output file which contains the selected entries and other useful information.
}


\examples{
library(Rindsel)
datos<-data.frame(read.csv("https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/C1_PSI_05_Phen.csv",header=T,na.strings=c(NA,".","-"))) #Raw data to analized.
datos$REP=as.factor(datos$REP)                                                                     #transform to factor.
datos$Block=as.factor(datos$Block)                                                                 #transform to factor.
datos$ENTRY=as.factor(datos$ENTRY)                                                                 #transform to factor.
file.wgt<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/weights_C1_PSI.csv"     #name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                    #Selection intensity.
design<-"lattice"                                                                            #Experimental design.
corr<-FALSE                                                                                  #You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
one.env<-TRUE                                                                                #Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                              #Use FALSE always.
softR<-""                                                                                    #Use "" always.
file.covG<-""                                                                                #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.
RESIMIndex(datos,file.wgt,selval,design,corr,out="outextRESIM.txt",outcsv="outRESIM.csv",rawdata,one.env,block.ex,softR,file.covG)
file.show("outextRESIM.txt")
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
