\name{LTIndex}
\Rdversion{1.0}
\alias{LTIndex}
\title{Linear Marker Selection Index}
\description{
Funtion for the Lande and Thompson selection index
}
\usage{
LTIndex(datos,file.wgt,selval,design,corr,out="outextLT.txt",outcsv="outLT.csv",
rawdata,file_nameMARK,file_nameQTL,one.env,block.ex,softR,file.covG)
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
  \item{file_nameMARK}{Text file in comma separated values (*.csv) format whose columns are markers and rows entries.}
  \item{file_nameQTL}{Text file in comma separated values (*.csv) format where we find QTL information (score and marker information for each trait).}
  \item{one.env}{Use FALSE for multienviromrent trials.}
  \item{block.ex}{Use FALSE always.}
  \item{softR}{Use "" always.}
  \item{file.covG}{When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.}
}

\details{
Lande and Thompson (1990) established the theoretical basis of the 
molecular marker-assisted selection (Molecular Marker Assisted 
Selection or MAS) for quantitative trait selection 
using computer simulation studies. Considerations 
Additional theoretical contributed more fully to 
understanding of fundamental aspects regarding the development of MAS 
population type, size n or sample and number of markers 
molecular (MM) required in MAS (Zhang and Smith 1992). 

The procedure of Lande and Thompson (1990) is based on a 
selection index which combines information 
with phenotypic information from the loci 
quantitative characters (Quantitative Trait Loci or QTL) associated with the MM. This is because 
it is not possible to identify all QTL affecting the nature of 
interest (Li, 1998), that is, unless all the QTLs affecting 
phenotypic variables of interest are identified, MAS 
phenotypic information should be combined with the effects 
of the QTL associated with MM, also known MQTL effects for 
ensure the efficiency of selection (Dekkers and Settar, 2004). In the 
construction of selection index Lande and Thompson (1990) 
requires: (1) identify the linkage between the GM and the QTL on a map 
of MMs, (2) estimate the effects MQTLs, and (3) combined with the effects MQTLs 
phenotypic information to classify individuals by 
a selection index, and subsequently develop 
lines, varieties or populations of interest. The effects MQTLs 
can be identified and estimated based on linkage disequilibrium 
created by crossing inbred lines or divergent populations 
(Jansen, 2003). 

As an illustrative example, assume that in the maize genome 
has a linkage map with ten molecular markers and the 
interest lies in enhancing the improvement of plant height; 
Suppose also that five QTLs with additive effects \eqn{ \alpha_1 
 } , \eqn{ \alpha_2 } ,{\ldots}, \eqn{ \alpha_5 }  have been identified in the genome 
associated with plant height and that the five QTLs are distributed 
at random from different chromosomes, ie, a QTL was found by 
chromosome. In practice, the idea of Lande and Thompson (1990) is 
to build a record (or score) by multiplying the values additives of the
QTLs for the coded values of molecular markers (MM) to 
they are bound, \eqn{ ie } , \eqn{ s = x_1 \alpha_1 + x_2 \alpha_2 + ... + x_5 \alpha_5 
 } , Where \eqn{x_i} , \eqn{i = 1,2, ..., } 5, denotes the coded values of the MM, 
which depend on the type of population that is used for 
selection in a population \eqn{F_2} the possible values \eqn{x_i} 
are 1, 0 and -1. In this situation, the IS is constructed as 
\deqn{ 
Y_{LT}=\beta_p p + \beta_s s 
} 
where \eqn{ p } denotes the height of plant, \eqn{s=x_1 \alpha_1 + x_2 \alpha_2 
+ ... + x_5 \alpha_5 }  and what is desired is to maximize the correlation between 
 \eqn{ Y_{LT} } and \eqn{ Z = \theta_p g_p } . In this case the vector can be shown that 
that maximizes the correlation between \eqn{Y_{LT} } y \eqn{Z=\theta _p g_p } is 
\eqn{beta_{LT} =T^{-1}K\theta_{LT} }, 
where \eqn{\beta '_{LT} =[\beta '& \beta'_s }, 
\eqn{\theta'_{LT} =[-1 & 0 ]}, \cr

T=                                      
\tabular{cc}{
\eqn{\sigma _p^2 }   \tab   \eqn{\sigma _s^2 }\cr
\eqn{\sigma _s^2 }   \tab   \eqn{\sigma _s^2 }\cr
}

K=                                      
\tabular{cc}{
\eqn{\sigma _g^2 }   \tab   \eqn{\sigma _s^2 }\cr
\eqn{\sigma _s^2 }   \tab   \eqn{\sigma _s^2 }\cr
}
where \eqn{\sigma _p^2 } is the phenotypic variance, 
 \eqn{ \sigma_g^2 } is the genotypic variance and \eqn{\sigma_s^2 } is the 
variance of \eqn{s=x_1 \alpha _1 +x_2 \alpha _2 +...+x_5 \alpha _5 }.

The case is considered more of a character is straightforward and 
All that changes is the size of the matrix \eqn{T} and 
\eqn{K}.
}

\value{
The LT selection index and an output file which contains the selected entries and other useful information.
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
file_nameMARK<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/C1_PSI_S2_05_Haplo.csv"  	#name of the file markers information.
file_nameQTL<-"https://raw.githubusercontent.com/RAngelaPG/RIndSel-R/main/data/QTL_scores_Simulated.csv"	#name of the file QTL information.
one.env<-TRUE                                                                                			 #Use FALSE for multienvironment trials.
block.ex<-FALSE                                                                              			 #Use FALSE always.
softR<-""                                                                                    			 #Use "" always.
file.covG<-""                                                                                			 #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.
LTIndex(datos,file.wgt,selval,design,corr,out="outextLT.txt",outcsv="outLT.csv",rawdata,file_nameMARK,file_nameQTL,one.env,block.ex,softR,file.covG)
file.show("outextLT.txt")
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
