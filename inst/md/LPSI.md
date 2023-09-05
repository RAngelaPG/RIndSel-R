#### Linear Phenotypic Selection Index LPSI (Smith).

In 1936, Smith developed the methodology for the selection simultaneously on several variables based on the phenotypic linear combination of

$Y = \beta'p$
$Z = \theta'g$

Where $Y$ is the selection index, $p = [p_{1} , p_{2},\ldots, p_{q}]$ and $\beta = [b_{1}, b_{2},\ldots, b_{q}]$ is the vector of values Phenotypic and the coefficient vector of $Y$ respectively, $Z$ denotes the breeding value or net genetic component that can be gained through various selection cycles, $g = [g_{1}, g_{2}, \ldots, g_{q}]$ is the vector of values genotypic, and $\theta = [\theta_{1}, \theta_{2},\ldots, \theta_{q}]$ is the vector of economic weights, which according to Smith (1936), can be determined according to the experience of researchers and is regarded as a vector of constants. With the aim establish the relationship between phenotypic values $P_{j}$ $(j = 1, 2,\ldots, q)$ observable and unobservable genotype in (Smith, 1936) proposed the following model 

$P_{j} = g_{j} + \epsilon_{j}$ 

Where $g_{j}$ is the genotypic value of the jth feature and $\epsilon_{j}$ is the environmental component, which affects the genotypic character. (Smith, 1936) assumes that the interactions between $g_{j }$ and $\varepsilon_{j}$ can be considered as a random effect, therefore $g_{j}$ has only additive effects such as $Z =\theta' g$ denotes the breeding value (Hazel, 1943; Kempthorne and Nordskog, 1959), under these assumptions, selection based on $Y= \beta' p$ leads to a response to selection (R) equal to 

$R=k \sigma_Z \rho_{YZ}=k \sigma_Z \frac{\theta' \Sigma \beta}{\sqrt{\theta'\Sigma\theta}\sqrt{\beta'S\beta}}$

Where $\Sigma$ and $S$ are matrices variance-covariance matrix for phenotypic and genotypic values,respectively, $k$ is the standard differential selection, $\theta'\Sigma \beta$ is the covariance between SI and the breeding value, $\sigma^2_Z=\theta'\Sigma\theta$ is the variance of $Z$ and $\sigma^2_Y=\beta'S \beta$ is the variance of SI, $\rho_{YZ}$ is the correlation between SI and the breeding value. In the Smith selection index $\beta_{s}=S^{-1} \Sigma \theta$ (where the subscript $S$ denotes the method of Smith, $S^{-1}$ is the inverse of the variance-covariance matrix of phenotypic $S$) allows us to construct the selection index, $Y =\beta_{s}'p$, which maximizes the correlation with the breeding value. 

```R

library(Rindsel)
datos<-read.csv("https://github.com/RAngelaPG/RIndSel-R/blob/master/data/C1_PSI_05_Phen.csv",header=T,na.strings=c(NA,"."."-")) #Raw data to analized.
file.wgt<-"https://github.com/RAngelaPG/RIndSel-R/blob/master/data/weigth_C1_PSI.csv")   #name of the file where we write the economic weights and restrictions. 
selval<-5                                                                                    #Selection intensity.
design<-"lattice"                                                                            #Experimental design.
corr<-FALSE                                                                                  #You can decide if you want to work with the correlation matrix instead of variance and covariance matrix.
rawdata<-TRUE                                                                                #By default is TRUE when you are using design option "lattice" or "rcbd", use FALSE for design option "AdjMeans".
one.env<-TRUE                                                                                #Use FALSE for multienviromrent trials.
block.ex<-FALSE                                                                              #Use FALSE always.
softR<-""                                                                                    #Use "" always.
file.covG<-""                                                                                #When design is "AdjMeans" and rawdata is FALSE, write the location of your variance and covariance matrix csv file.

LPSI(datos,file.wgt,selval,design,corr,out="outextBLPSI.txt",outcsv="outBLPSI.csv",rawdata,one.env,block.ex,softR,file.covG)

```
