\name{SmithIndex}
\Rdversion{1.0}
\alias{SmithIndex}
\title{Linear Phenotypic Selection Index}
\description{
Funtion for Smith Selection Index
}
\usage{
SmithIndex(datos,file.wgt,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",
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
In 1936, Smith developed the methodology for the selection simultaneously on several variables based on the phenotypic linear combination of

\deqn{Y = \beta'p} 
\deqn{Z = \theta'g} 

Where \bold{Y} is the selection index (SI), \eqn{p = [p_1, p_2, \ldots, p_q]} and \eqn{\beta = [b_1, b_2, \ldots, b_q]} is the vector of values phenotypic and the coefficient vector of \bold{Y} respectively, \bold{Z} denotes the breeding value or net genetic component that can be gained through various selection cycles, \eqn{g = [g_1, g_2,\ldots, g_q]} is the vector of values genotypic, and \eqn{\Theta = [\Theta_1, \Theta_2, \ldots, \Theta__q]} is the vector of economic weights, which according to Smith (1936), can be determined according to the experience of researchers and is regarded as a vector of constants. With the aim 
establish the relationship between phenotypic values \eqn{P_j} (j = 1, 2,..., q) observable and unobservable genotype in (Smith, 1936) proposed the following model 
\deqn{P_j = g_j + \epsilon_j} 

Where \bold{g_j} is the genotypic value of the jth feature and \eqn{\epsilon_j} is the environmental component, which affects the genotypic character. (Smith, 1936) assumes that the interactions between \bold{g_j} and \eqn{\varepsilon_j} can be considered as a 
random effect, therefore \bold{g_j} has only additive effects such as \eqn{Z = \Theta'g} denotes the breeding value (Hazel, 1943; Kempthorne and Nordskog, 1959), under these assumptions, selection based on \eqn{Y = \beta'p} leads to a response to selection (R) equal to 
\deqn{R = k \sigma_Z \rho_YZ = k \sigma_Z (\theta' \Sigma \beta) / (\sqrt \theta'\Sigma\theta \sqrt \beta'S\beta )}

Where \bold{\eqn{\Sigma}} and \bold{S} are matrices variance-covariance matrix for phenotypic and genotypic values, respectively, \bold{k} is the standard differential selection, \bold{ \eqn{ \theta } ' \eqn{ \Sigma \beta } } is the covariance between \bold{Y}, \bold{Z}, \bold{\eqn{\theta }'\eqn{\Sigma\theta}} is the variance of \bold{Z} and \bold{ \eqn{ \beta } 'S \eqn{ \beta } } is the variance of \bold{Y}, \eqn{\rho_yz} is the correlation between \bold{Y} and \bold{Z}. 
In the Smith selection index \eqn{\beta_s = S^-1 \Sigma\theta} (where the subscript \bold{S} denotes the method of Smith, \eqn{S^-1} is the inverse of the variance-covariance matrix of phenotypic \bold{S}) allows us to construct the selection index (SI), \eqn{Y = \beta_s'p}, which maximizes the correlation with the breeding value. \cr

}

\value{
The Smith's selection index and an output file which contains the selected entries and other useful information.
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
corr<-FALSE                                                                                        #You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                      #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
one.env<-TRUE                                                                                      #Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                                    #Use FALSE always.
softR<-""                                                                                          #Use "" always.
file.covG<-""                                                                                      #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.
SmithIndex(datos,file.wgt,selval,design,corr,out="outextLPSI.txt",outcsv="outLPSI.csv",rawdata,one.env,block.ex,softR,file.covG)
file.show("outextLPSI.txt")
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
