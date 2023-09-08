\name{MESIMIndex}
\Rdversion{1.0}
\alias{MESIMIndex}
\title{Molecular Eigen Selection Index Method}
\description{
Funtion for MESIM selection index
}
\usage{
MESIMIndex(datos,file.wgt,selval,design,corr,method,out="outextLT.txt",outcsv="outLT.csv",
rawdata,file_nameMARK,one.env,block.ex,softR,file.covG)
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
  \item{file_nameMARK}{Text file in comma separated values (*.csv) format whose columns are markers and rows entries.}
  \item{one.env}{Use FALSE for multienviromrent trials.}
  \item{block.ex}{Use FALSE always.}
  \item{softR}{Use "" always.}
  \item{file.covG}{When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.}}

\details{
MESIM is a generalization of ESIM to the case in which is incorporated 
MM information to IS, similarly as in the IS Lande and 
Thompson (1990). Following the basic idea of Kempthorne and Nordskog 
(1959), Ceron-Rojas {et al}. (2008b), maximized the correlation between \eqn{ Y_M 
={\rm{\bf{\beta} '}}_p{\rm{\bf p}} +{\rm{\bf{\beta}'}}_s{\rm{\bf 
s}} ={\rm{\bf{\beta}'}}_M{\rm{\bf p}}_{ps} }{ Y_M =\beta'_p*P +\beta'_s*S =\beta'_M*P_{ps} } and \eqn{Z_M ={\rm{\bf{\theta 
} '}}_1{\rm{\bf g}} +{\rm{\bf{\theta}'}}_2{\rm{\bf s}} ={\rm{\bf 
{\theta}'}}_M{\rm{\bf g}}_{gs}}{Z_M =\theta'_1*g +\theta_2 *s=\theta'_M*g_{gs}} , \eqn{\rho_{Y_M Z_M}^2} , where \eqn{{\rm{\bf 
p}} }{p} and \eqn{ {\rm{\bf g}} }{g} have been defined in Equations 1 and 2, and \eqn{ {\rm 
{\bf s}} }{s} is the vector of records of the additive effects of QTLs 
associated with MM, \eqn{ {\rm{\bf{p} '}}_{ps}=\left [{{\begin{array}{*{20} c} 
{{\rm{\bf{p} '}}} \hfill &{{\rm{\bf{s}'}}} \hfill \\
\end{array}}} \right] }{P'_{ps}=[p'  s']} , \eqn{ {\rm{\bf{g} '}}_{gs}=\left [ 
{{\begin{array}{*{20} c} 
{{\rm{\bf{g} '}}} \hfill &{{\rm{\bf{s}'}}} \hfill \\
\end{array}}} \right] }{G'_{gs}=[g'  s']} , \eqn{ {\rm{\bf{\beta} '}}_M = \left [ 
{{\begin{array}{*{20} c} 
{{\rm{\bf{\beta} '}}_P} \hfill &{{\rm{\bf{\beta}'}}_s} \hfill \\
\end{array}}} \right] }{\beta'_M=[\beta'_P   \beta'_s]} , and \eqn{ {\rm{\bf{\theta} '}}_M = \left [ 
{{\begin{array}{*{20} c} 
{{\rm{\bf{\theta} '}}_1} \hfill &{{\rm{\bf{\theta}'}}_2} \hfill \\
\end{array}}} \right] }{\theta'_M=[\theta'_1   \theta'_2]} . In the illustrative example of Lande and 
Thompson (1990), \eqn{ Y_M = \beta_p*p + \beta_s*s } , and \eqn{ Z_M = \theta_p*G_P + \theta 
_s*s } , where \eqn{ p } denotes the plant height and \eqn{ s = x_1 \alpha_1 + x_2 \alpha_2 + ... + x_5 \alpha_5}. 

Again, because \eqn{ \rho_{Y_M Z_M }^2}  is invariant to changes in 
is maximized scale \eqn{\rho_{Y_M Z_M ^2} } under the constraints 
 \eqn{ {\rm{\bf{\beta} '}}_M{\rm{\bf T \beta}}_M = 1 }{ \beta '_M*T*\beta_M = 1 } and \eqn{ {\rm{\bf{\theta 
} '}}_M{\rm{\bf K \theta}}_M =1 }{ \theta '_M*K*\theta_M = 1 }, so in MESIM, it is necessary 
maximize 
\deqn{ 
\Phi = \left ({{\rm{\bf{\theta} '}}_M{\rm{\bf K \beta}}_M} \right) ^ 2 - \mu 
\left ({{\rm{\bf{\beta} '}}_M{\rm{\bf T \beta}}_M -1} \right) - \omega 
\left ({{\rm{\bf{\theta} '}}_M{\rm{\bf K \theta}_M} -1} \right) 
}{ (\theta '_M*K*\beta_M-1 )^2-\mu(\beta '_M*T*\beta_M-1) -\omega(\theta '_M*K*\theta_M-1 )} 
With respect to \eqn{ {\rm{\bf \beta}}_M }{\beta_M} , \eqn{ {\rm{\bf \theta}}_M }{\theta_M} , \eqn{ \mu } , and 
 \eqn{ \omega } , where \eqn{ {\rm{\bf \beta}}_M }{\beta_M} is the vector of coefficients 
 \eqn{ \mbox{MESIM} }{MESIM} , \eqn{ {\rm{\bf \theta}}_M }{\theta_M} is the vector of coefficients 
 \eqn{ Z_M ={\rm{\bf{\theta} '}}_M{\rm{\bf g}}_{gs} }{ Z_M =\theta'_M* g_{gs}} and \eqn{ \mu } and \eqn{ \omega } 
are Lagrange multipliers. Ceron-Rojas {et al}. (2008b) found that 
the solution is equal 
\deqn{ 
\label{EQ4} 
({\rm{\bf Q}} - \mu{\rm{\bf I}}){\rm{\bf \beta}}_M ={\rm{\bf 0}}, 
}{(Q-\mu I)\beta_M=0} 
where \eqn{ {\rm{\bf Q}}={\rm{\bf T}} ^{-1}{\rm{\bf K}} }{Q=T^{-1}*K} . Thus, in 
 \eqn{ \mbox{MESIM} }{}MESIM , the value that maximizes \eqn{ \rho_{Y_M Z_M} ^ 2 } under the 
constraints \eqn{ {\rm{\bf{\beta} '}}_M{\rm{\bf T \beta}}_M = 1 }{\beta '_M*T*\beta_M=1} and \eqn{ {\rm 
{\bf{\theta} '}}_M{\rm{\bf K \theta}}_M = 1}{\theta '_M*K*\theta_M=1} is the first eigenvalue ( \eqn{ \mu 
) } Matrix \eqn{ {\rm{\bf Q}} }{Q} , and the vector for building  \eqn{ Y_M } 
(With maximum correlation with \eqn{ Z_M ={\rm{\bf{\theta}'}}_M{\rm{\bf 
g}}_{gs}) }{ Z_M =\theta'_M* g_{gs}} is the first eigenvector ( \eqn{ {\rm{\bf \beta}}_M) }{\beta_M}
matrix 
 \eqn{ {\rm{\bf Q}} }{Q} .
}

\value{
The MESIM selection index and an output file which contains the selected entries and other useful information.
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
method<-"vanraden"																				   #For create the relationship marker matrix.
rawdata<-TRUE                                                                                      #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
file_nameMARK<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/C1_PSI_S2_05_Haplo.csv")   #name of the file markers information.
one.env<-TRUE                                                                                      #Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                                    #Use FALSE always.
softR<-""                                                                                          #Use "" always.
file.covG<-""                                                                                      #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.
MESIMIndex(datos,file.wgt,selval,design,corr,method,out="outextMESIM.txt",outcsv="outMESIM.csv",rawdata,file_nameMARK,one.env,block.ex,softR,file.covG)
file.show("outextMESIM.txt")
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
