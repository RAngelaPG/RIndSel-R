\name{ESIMIndex}
\Rdversion{1.0}
\alias{ESIMIndex}
\title{Eigen Selection Index Method.}
\description{
This function estimates the selection index by the Eiegen Selection Index Method (ESIM).
}
\usage{
ESIMIndex(datos,file.wgt,markers=NULL,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",
rawdata, one.env, block.ex, softR, file.covG)
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
The eigen selection index method (or MESIM for its acronym in English) is a selection method based on the IS Smith (1936) in which
uses the theory of singular value decomposition (DS) to estimate the vector \eqn{\beta} that maximizes the correlation \eqn{\rho_ZY}. In this case the first eigenvector elements (\eqn{\beta}) of \eqn{P^-1G} determine the proportion to each character contributes to IS \eqn{Y = \beta'p} , also the first eigenvalue of \eqn{P^-1G} is used in the response to selection, such IS was proposed by Ceron-Rojas {et al}. (2008a, 2008b). 

MESIM theory is directly related to the theory of the canonical correlations in which, according to Muirhead (2005), when a set of variables may be divided so natural two subsets of random variables, and what you want is maximize the correlation between the two subsets of variables, canonical correlation analysis is very efficient. Consider the set of phenotypic variables (\eqn{p}) and the set of genotypic variables (\eqn{g}) in the context the improvement of plant and animal, in this case must be, indeed, the vector of phenotypic variables, \eqn{p} , and the vector of genotypic variables, \eqn{g} , belong, naturally, two sets of variables. Analysis canonical correlations reduces the correlation between \eqn{p} and \eqn{g} to its simplest form by linear transformations of \eqn{p} and \eqn{g} , ie  \eqn{Y = \beta'p} and \eqn{Z = \theta'g} . 

Using a concept similar to that of Kempthorne and Nordskog (1959), Ceron-Rojas {et al}. (2008a, 2008b) maximized the response to selection (Equation 2) to maximize \eqn{\rho_YZ^2} . Note that the variances  \eqn{Y = \beta'p} and \eqn{Z = \theta'g} are constant each selection cycle, therefore, the selection of genotypes can be done using \eqn{Y = \beta'p} or \eqn{Y / \sqrt \beta'P\beta} . Because of this, by maximizing \eqn{\rho_YZ^2} is possible to introduce the constraints \eqn{\beta'P\beta = 1} and \eqn{\theta'G \theta = 1}  so that in MESIM is necessary to maximize \eqn{\Phi = (\theta'G \beta)^2 - \mu  (\beta'P \beta -1 ) - \omega  (\theta'G \theta -1 )} with respect to \eqn{\beta} , \eqn{\theta} , \eqn{ \mu } , and  \eqn{ \omega } , where \eqn{\beta} is the vector of coefficients  \eqn{ Y = \beta'p } , \eqn{\theta} is the vector coefficients of \eqn{Z = \theta'g} , and \eqn{ \mu } and \eqn{ \omega } are multipliers Lagrange. In MESIM the values of \eqn{\theta } does not necessarily economic weights. 

When \eqn{ \Phi } is derived with respect to \eqn{\beta} and \eqn{\theta }, and the result is equal to zero vector, we have 
\deqn{(\theta 'G \beta) ^ 2 G \theta - \mu P \beta = 0} 
\deqn{(\theta 'G \beta) ^ 2 G \beta - \omega G \theta =0}
 
By the constraints \eqn{\beta'P \beta = 1 } and \eqn{\theta'G \theta = 1} ,when the last equation  is multiplied by \eqn{\beta'} and the next equation is multiplied by \eqn{\theta'} ,the result is \eqn{({\theta}'G \beta)^2 = \omega = \mu }. Therefore, \eqn{ \mu } maximizes \eqn{\rho_YZ^2} under the constraints \eqn{\beta'P \beta =1 } and \eqn{\theta'G \theta = 1}. 

The next problem is to determine the vector \eqn{\beta } that allows building the IS \eqn{Y =\beta'p} that has maximum correlation with \eqn{Z=\theta'g } . According to Ceron-Rojas {et al}. (2008a, 2008b), the \eqn{ \beta } required is the solution to the  following equation 
\deqn{( Q - \mu I)\beta = 0,} where \eqn{ Q= P^-1 G } . 

Thus, in MESIM , the value that maximizes \eqn{\rho_{YZ}^2 } under the constraints \eqn{ \beta'P \beta=1 } and \eqn{ \theta'G \theta=1 } is the first eigenvalue ( \eqn{ \mu) } matrix \eqn{ Q } , and the vector allows to construct \eqn{Y=\beta'p } (with maximum correlation with \eqn{ Z=\theta'g }) is the first eigenvector ( \eqn{ \beta})  Matrix \eqn{Q}.

}

\value{
The ESIM selection index and an output file which contains the selected entries and other useful information.
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


\examples{
library(Rindsel)
datos<-data.frame(read.csv("https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/C1_PSI_05_Phen.csv",header=T,na.strings=c(NA,".","-"))) #Raw data to analized.
datos$REP=as.factor(datos$REP)																		                                #transform to factor.
datos$Block=as.factor(datos$Block)																	                              #transform to factor.
datos$ENTRY=as.factor(datos$ENTRY)																	                              #transform to factor.
file.wgt<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/weights_C1_PSI.csv"   #name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                        #Selection intensity.
design<-"lattice"                                                                                #Experimental design.
corr<-FALSE                                                                                      #You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                    #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
one.env<-TRUE                                                                                    #Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                                  #Use FALSE always.
softR<-""                                                                                        #Use "" always.
file.covG<-""                                                                                    #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.
ESIMIndex(datos,file.wgt,markers=NULL,selval,design,corr,out="outextESIM.txt",outcsv="outESIM.csv",rawdata,one.env,block.ex,softR,file.covG)
file.show("outextESIM.txt")
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

\keyword{ package }
