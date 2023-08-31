MESIMIndex <-function(file.dat=NULL,file.wgt=NULL,selval=5,design="lattice",corr=FALSE,method=c("vanraden","classic"),
				out="out.txt",outcsv="out.csv",rawdata=TRUE,markers=NULL,one.env=T, block.ex, softR, file.covG){
    if(rawdata==TRUE){
    traits <<- NULL
    traits2 <<- NULL
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
    traits <<- traits <<- Read(file.dat, file.wgt,design,mis =  c(".","NA","na","N/A"), one.env, genomic=TRUE, markers=markers, block.ex, softR)
	nom.traits <- colnames(traits)[-(1:fcol - 1)]
  }
  
  
  theta <- rep(wgts[match(nom.traits,wgts[,1]),2],2)
  names(theta)<-rep(nom.traits,2)
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
  #QTLS <- ReadQTL()
  cat("\n","   Wait ...","\n\n")
  cov.scores <- GEBV(P=covP,G=covG,X=codes,method=method)
  MVMNS <<- cov.scores$ga 
  dimnames(MVMNS) <- dimnames(covG)
  
  if(corr==TRUE) {
    mat.name <- "CORRELATION"
    MVGS <<- cov2cor(covG)
    MVG <- MVGS
    MVPS <<- cov2cor(covP)
    MVP <- MVPS
    MVMS <<- cov2cor(MVMNS)
    MVM <- MVMS
    
  } else{
    mat.name <- "COVARIANCE"
    MVG <- covG
    MVP <- covP
    MVM <- MVMNS
  }
  
 
  MMVP <<- rbind(cbind(MVP,MVM),cbind(MVM,MVM))
  MMVG <<- rbind(cbind(MVG,MVM),cbind(MVM,MVM))
  svdIPG <- svd(solve(MMVP)%*%MMVG)
  EigVecMESIM <- as.matrix(abs(svdIPG$u[,1]))*theta
  EigVecMESIM <- EigVecMESIM/as.vector(sqrt(t(EigVecMESIM)%*%EigVecMESIM))
  scores <- cov.scores$ebv
  YMESIM <<- cbind(scale(allmean),scale(scores))%*%EigVecMESIM
  rnames <- rownames(YMESIM)
  thetaMESIM <- solve(MMVG)%*%(MMVP)%*%EigVecMESIM
  thetaNMESIM <- thetaMESIM%*%(1/sqrt(t(thetaMESIM)%*%thetaMESIM))
  vGb <- t(thetaNMESIM)%*%MMVG%*%EigVecMESIM
  VIDS <- t(EigVecMESIM)%*%MMVP%*%EigVecMESIM
  VDCG <- t(thetaNMESIM)%*%MMVG%*%thetaNMESIM
  vGv <- sqrt(abs(VDCG))
  bPb <- sqrt(abs(VIDS))
  corrAnalysis <- min(0.9999,vGb/(vGv%*%bPb))  
  selected <- YMESIM>=quantile(YMESIM,1-selval/100)
  selYMESIM <- YMESIM[selected]
  ord <- order(selYMESIM,decreasing=TRUE)
  selYMESIM <- selYMESIM[ord]
  selentry <-  allmean[selected,]
  nom.selentry <- rownames(allmean)[selected]
  selentry <- data.frame(as.matrix(selentry)[ord,])
  dimnames(selentry) <- list(nom.selentry[ord],nom.traits)
  MMESIMSI <- apply(selentry,2,mean)
  MMESIMall <- apply(allmean,2,mean)
  GAINbyMEANS=MMESIMSI-MMESIMall
  ks <-(100/selval)*(1/sqrt(2*pi))*exp((-qnorm(1-(selval/100))^2)/2)
  Rsel<-PlotResp(t(EigVecMESIM)%*%MMVP%*%EigVecMESIM,"MESIM",selval)
  #Rsel<-ks*sqrt(t(EigVecMESIM)%*%MMVP%*%EigVecMESIM)
  H2<-(t(EigVecMESIM)%*%MMVG%*%EigVecMESIM)/(t(EigVecMESIM)%*%MMVP%*%EigVecMESIM)
  MESIMGain <- as.vector(ks*(MMVG%*%EigVecMESIM)/as.numeric(sqrt(t(EigVecMESIM)%*%MMVP%*%EigVecMESIM)))
  selentry2 <- rbind(selentry,MMESIMSI,MMESIMall,GAINbyMEANS,MESIMGain)
  selMESIM <- data.frame(selentry2,c(selYMESIM,rep(NA,4)))
  colnames(selMESIM)[ncol(selMESIM)] <- "MESIM index"
  rnames <- c("Mean of Selected Individuals", "Mean of all Individuals",
            'Selection Differential',paste("Expected Genetic Gain for",paste(selval,"%",sep=""),sep=" "))
  rownames(selMESIM)[(length(selYMESIM)+1):(length(selYMESIM)+4)] <- rnames
  allMESIM <- data.frame(allmean,YMESIM)
   colnames(allMESIM)[ncol(allMESIM)] <- "MESIM index"
  cat("MESIM SELECTION INDEX METHOD","\n",file=out)
  cat("\n",paste("GENETIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(MVG,2),file=out,col.names=T,append=T)
  cat("\n",paste("PHENOTYPIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(MVP,2),file=out,col.names=T,append=T)
  cat("\n","MOLECULAR COVARIANCE MATRIX","\n",file=out,append=T)
  print.char.matrix(round(MVM,2),file=out,col.names=T,append=T)
  if (corr != TRUE) {
  cat("\n\n","COVARIANCE BETWEEN THE MESIM SELECTION INDEX AND THE BREEDING VALUE: ", 
      round(vGb,3),"\n",file=out,append=T)
  cat("\n","VARIANCE OF THE MESIM SELECTION INDEX:                                 ", 
      round(VIDS,3),"\n",file=out,append=T)
  cat("\n","VARIANCE OF THE BREEDING VALUE:                                        ", 
      round(VDCG,3),"\n",file=out,append=T)
  cat("\n","CORRELATION BETWEEN THE MESIM SELECTION INDEX AND THE BREEDING VALUE:  ", 
      round(corrAnalysis,3),"\n",file=out,append=T)
  cat("\n", "RESPONSE TO SELECTION:                                                ", 
      round(Rsel,3), "\n", file = out, append = T)
  cat("\n", "HERITABILITY:                                                         ", 
      round(H2,3), "\n", file = out, append = T)
  }	  
  cat("\n\n",
      paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE MESIM SELECTION, MEANS AND GAINS FOR", paste(selval,"%",
            sep =""),sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(selMESIM,2),file=out,col.names=T,append=T)
  cat("\n\n","VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE MESIM SELECTION INDEX","\n",file=out,append=T)
  print.char.matrix(round(allMESIM,2),file=out,col.names=T,append=T)
  #file.show(out)
  write.csv(selentry,na="",file=outcsv)
  #write.csv(allMESIM,file=paste("all",outcsv,sep=""))
  allprint=allMESIM
  if(rawdata == TRUE){SummaryStats(traits, design, softR, one.env, block.ex, covG, covP,allprint, "BLUE_BLUP_H2", nom.traits)}
  save(MVG,MVP,vGb,VIDS,VDCG,corrAnalysis,Rsel,H2,MVM,selentry,allprint,file="forReport.RData")
  #cat("\n",paste("Output saved in the files:",paste(out,",",sep=""),outcsv,"and",paste("all",outcsv,sep=""),sep=" "),"\n")
  #rm(fcol, ntraits, mat, pos = ".GlobalEnv")  
}

