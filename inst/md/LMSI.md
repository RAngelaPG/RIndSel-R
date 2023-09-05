## Linear Marker Selection Index (LMSI-LT).

Lande and Thompson (1990) established the theoretical basis of the molecular marker-assisted selection (Molecular Marker Assisted Selection or MAS) for quantitative trait selection using computer simulation studies. Considerations additional theoretical contributed more fully to understanding of fundamental aspects regarding the development of MAS population type, size \~{n} or sample and number of markers molecular (MM) required in MAS (Zhang and Smith 1992). 

The procedure of Lande and Thompson (1990) is based on a selection index which combines information with phenotypic information from the loci quantitative characters (Quantitative Trait Loci, or QTL) associated with the MM. This is because it is not possible to identify all QTL affecting the nature of interest (Li, 1998), that is, unless all the QTLs affecting phenotypic variables of interest are identified, MAS phenotypic information should be combined with the effects of the QTL associated with MM, also known MQTL effects for ensure the efficiency of selection (Dekkers and Settar, 2004). In the construction of selection index Lande and Thompson (1990) requires: (1) identify the linkage between the GM and the QTL on a map of MMs, (2) estimate the effects MQTLs, and (3) combined with the effects MQTLs phenotypic information to classify individuals by a selection index, and subsequently develop lines, varieties or populations of interest. The effects MQTLs can be identified and estimated based on linkage disequilibrium created by crossing inbred lines or divergent populations (Jansen, 2003). 

As an illustrative example, assume that in the maize genome has a linkage map with ten molecular markers and the interest lies in enhancing the improvement of plant height; Suppose also that five QTLs with additive effects $\alpha_1$ , $\alpha_2$ ,$\ldots$, $\alpha_5$ have been identified in the genome associated with plant height and that the five QTLs are distributed at random from different chromosomes, ie, a QTL was found by chromosome. In practice, the idea of Lande and Thompson (1990) is to build a record (or score) by multiplying the values additives of the QTLs for the coded values of molecular markers (MM) to they are bound, ie, $s = x_1 \alpha_1 + x_2 \alpha_2 + ... + x_5 \alpha_5$ , Where $x_i$ , $i = 1,2, ...,5$, denotes the coded values of the MM, which depend on the type of population that is used for selection in a population $F_2$ the possible values $x_i$ are 1, 0 and -1. In this situation, the IS is constructed as $Y_{LT}=\beta_p p + \beta_s s$ where $p$ denotes the height of plant, $s=x_1 \alpha_1 + x_2 \alpha_2+ ... + x_5 \alpha_5$ and what is desired is to maximize the correlation between $Y_{LT}$ and $Z = \theta_p g_p$ . In this case the vector can be shown that maximizes the correlation between $Y_{LT}$ y $Z=\theta _p g_p$ is $\beta_{LT} =T^{-1}K\theta_{LT}$, where $\beta '_{LT} =(\beta \beta'_s)$, $\theta'_{LT} =(-1 & 0)$, $T=$ ```math \begin{bmatrix}$\sigma _p^2$ & $\sigma _s^2$\\$\sigma _s^2$ &  $\sigma _s^2$\end{bmatrix}``` and  $K=$ ```math \begin{bmatrix}$\sigma _g^2$ & $\sigma _s^2$\\$\sigma _s^2$ & $\sigma _s^2$\end{bmatrix}``` where $\sigma _p^2$ is the phenotypic variance, $\sigma_g^2$ is the genotypic variance and $\sigma_s^2$ is the variance of $s=x_1 \alpha _1 +x_2 \alpha _2 +...+x_5 \alpha _5$.

The case is considered more of a character is straightforward and all that changes is the size of the matrix $T$ and $K$.

```R

library(Rindsel)
datos<-read.csv("https://github.com/RAngelaPG/RIndSel-R/blob/master/data/C1_PSI_05_Phen.csv",header=T,na.strings=c(NA,"."."-")) #Raw data to analized.
file.wgt<-"https://github.com/RAngelaPG/RIndSel-R/blob/master/data/weigth_C1_PSI.csv")   #name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                    #Selection intensity.
design<-"lattice"                                                                            #Experimental design.
corr<-FALSE                                                                                  #You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
file_nameMARK<-
file_nameQTL<-
one.env<-TRUE                                                                                #Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                              #Use FALSE always.
softR<-""                                                                                    #Use "" always.
file.covG<-""                                                                                #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.

LTIndex(datos,file.wgt,selval,design,corr,out="outextLT.txt",outcsv="outLT.csv",rawdata,file_nameMARK,file_nameQTL,one.env,block.ex,softR,file.covG)

```
[Return to examples](https://github.com/RAngelaPG/RIndSel-R/blob/master/Readme.md)
