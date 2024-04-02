LGSI<-function(file.dat=NULL,file.wgt=NULL,selval=5,design="lattice",corr=FALSE,method="vanraden",out="out.txt",outcsv="out.csv",rawdata=TRUE,markers=NULL,one.env=T,block.ex, softR, file.covG,Tmarkers=NULL,LG=1)
{
  if(rawdata==TRUE){
    traits <<- NULL
    covG <<- NULL
    covP <<- NULL
    allmean <<- NULL
  }else {
      traits <<- file.dat
      nom.traits <- colnames(traits)[-1]
	  wgts <<- read.csv(file = file.wgt, header = T)  
	  wgts <<-as.data.frame(wgts)
      covG <<- as.matrix(read.csv(file = file.covG, header = T)[,-1])
      dimnames(covG) <- list(nom.traits, nom.traits)
      allmean <<- traits[ ,-1]
      colnames(allmean) <- nom.traits
      rownames(allmean) <- file.dat[,1]
      covP <<- fCovariances(allmean)
      dimnames(covP) <- dimnames(covG)
    }
	
  if(is.null(traits)){
    traits <<- Read(file.dat, file.wgt,design,mis =  c(".","NA","na","N/A"), one.env, genomic=TRUE, markers=markers, block.ex, softR)
	nom.traits <- colnames(traits)[-(1:fcol - 1)]
  }
    
  theta <- wgts[match(nom.traits,wgts[,1]),2]
  names(theta)<-nom.traits
  n<-length(theta)
  cat("\n","   Wait ...","\n\n")  
  if(is.null(covG)&is.null(covP)&is.null(allmean)){
    covG <<- gCovariances(traits, design, softR,one.env,block.ex) 
    dimnames(covG) <- list(nom.traits, nom.traits)
    allmean <<- Means(traits,design, softR,one.env,block.ex)
    colnames(allmean) <- nom.traits
    covP <<- fCovariances(allmean)
    dimnames(covP) <- dimnames(covG)
  }
  
  codes <- ReadMarkers(markers)  
  cat("\n","   Wait ...","\n\n")
  cov.scores <- GEBVLGSI(P=covP,G=covG,X=codes,method=method)
  MVMNS <<- cov.scores$ga 
  dimnames(MVMNS) <- dimnames(covG)
  
  Wt <<- ReadTestingLGSI(Tmarkers)  
  gebv<-Wt[,-1]%*%cov.scores$ebvs 
  varGamma<-1/nrow(Wt[,-1])*t(scale(gebv, scale=FALSE))%*%(1/cov.scores$cm*(Wt[,-1]%*%t(Wt[,-1])))%*%scale(gebv, scale=FALSE)
  #varmol<-varMol<-cov.scores$Vga
  if(corr==TRUE){
    mat.name <- "CORRELATION"
    MVGS <<- cov2cor(covG)
    MVG <- MVGS
    MVPS <<- cov2cor(covP)
    MVP <- MVPS
    MVMS <<- cov2cor(MVMNS)
    MVM <- varGamma
	dimnames(MVM) <- list(nom.traits,nom.traits) 
  }else{
    mat.name <- "COVARIANCE"
    MVG <- covG
    MVP <- covP
    MVM <- varGamma
	dimnames(MVM) <- list(nom.traits,nom.traits) 
  }
  
  MMVP <<- MVM 
  MMVG <<- MVG 
  svdIPG <- svd(solve(MMVG)%*%MMVP)   
  #LG<-1  
  Yb <<- gebv%*%theta  
  rownames(Yb)<-paste("Entry",as.character(Wt[,1]))
  colnames(Yb)<-"LGSI"    
  vGb <- t(theta)%*%MMVG%*%theta
  VIDS <- t(theta)%*%MMVP%*%theta
  VDCG <- t(theta)%*%MMVG%*%theta
  vGv <- sqrt(abs(VDCG))
  bPb <- sqrt(abs(VIDS))
  corrAnalysis <- min(0.9999,vGb/(vGv%*%bPb))  
  selected <- Yb>=quantile(Yb,1-selval/100)
  ord <- order(Yb[selected],decreasing=TRUE)
  selentry <-  Yb[selected,1]
  selentry <- data.frame(as.matrix(selentry)[ord,])
  colnames(selentry)<-"LGSI"
  ks <-(100/selval)*(1/sqrt(2*pi))*exp((-qnorm(1-(selval/100))^2)/2)  #selection intensity
  Rsel1<-(ks/LG)*(vGb/sqrt(t(theta)%*%MMVP%*%theta))
  Rsel<-PlotResp(Rsel1/(ks/LG),"LGSI",selval)
  H2<-t(theta)%*%MMVG%*%theta/t(theta)%*%MMVP%*%theta
  Gain <- (1/sqrt(t(theta)%*%MMVG%*%theta))*as.vector((ks/LG)*theta%*%MMVP)
  Gain <- as.data.frame(t(Gain))
  names(Gain)<-colnames(MMVG)
  cat("LGSI SELECTION INDEX METHOD","\n",file=out)
  cat("\n",paste("GENETIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(MVG,2),file=out,col.names=T,append=T)
  cat("\n",paste("PHENOTYPIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(MVP,2),file=out,col.names=T,append=T)
  cat("\n","MOLECULAR COVARIANCE MATRIX","\n",file=out,append=T)
  print.char.matrix(round(MVM,2),file=out,col.names=T,append=T)  
  if (corr != TRUE) {
  cat("\n\n","COVARIANCE BETWEEN THE LGSI SELECTION INDEX AND THE BREEDING VALUE: ", round(vGb,3),"\n",
      file=out,append=T)
  cat("\n","VARIANCE OF THE LGSI SELECTION INDEX:                                 ", round(VIDS,3),"\n",
      file=out,append=T)
  cat("\n","VARIANCE OF THE BREEDING VALUE:                                       ",round(VDCG,3),"\n",
      file=out,append=T)
  cat("\n","CORRELATION BETWEEN THE LGSI SELECTION INDEX AND THE BREEDING VALUE:  ",round(corrAnalysis,3),"\n",
      file=out,append=T)
  cat("\n", "RESPONSE TO SELECTION FOR TESTING POPULATION:                        ", 
      round(Rsel1,3), "\n", file = out, append = T)
  cat("\n", "HERITABILITY FOR TESTING POPULATION:                                 ", 
      round(H2,3), "\n", file = out, append = T)
  cat("\n", paste("EXPECTED GENETIC GAIN PER TRAIT", sep = " "),"\n", file = out, append = T)
  print.char.matrix(round(Gain, 3), file = out, col.names = T, append = T) 
   } 
  cat("\n\n",paste("VALUE OF THE LGSI SELECTION", paste(selval,"%",sep =""),sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(selentry,2),file=out, col.names=T,append=T)
  cat("\n\n","VALUES OF THE LGSI SELECTION INDEX",
      "\n",file=out,append=T)
  print.char.matrix(round(Yb,2),file=out,col.names=T,append=T)
allprint=Yb
if(rawdata == TRUE){SummaryStats(traits, design, softR, one.env, block.ex, covG, covP,allprint, "BLUE_BLUP_H2", nom.traits)}
  save(MVG,MVP,vGb,VIDS,VDCG,corrAnalysis,Rsel,H2,MVM,selentry,allprint,file="forReport.RData")
  #file.show(out)  
}
