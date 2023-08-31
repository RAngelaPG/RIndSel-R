LTIndex <-function(file.dat=NULL,file.wgt=NULL,selval=5,design="lattice",
    corr=FALSE,out="out.txt",outcsv="out.csv",rawdata=TRUE,markers=NULL,qtl=NULL,one.env=T, block.ex, softR, file.covG){
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
  
  
  theta <- rep(wgts[match(nom.traits,wgts[,1]),2],2)
  theta1 <- rep(wgts[match(nom.traits,wgts[,1]),2],1)
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
  QTLS <- ReadQTL(qtl)
  cov.scores <- CovMol(codes,QTLS)
  MVMNS <<- as.matrix(cov.scores$S[1:ncol(allmean),1:ncol(allmean)])
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
  scores <- cov.scores$Scores[,1:ncol(allmean)]
  BLT1 <<- as.vector(solve(MMVP)%*%MMVG%*%theta)
  BLT2 <<- as.vector(sqrt(t(BLT1)%*%BLT1))
  BLT <<- BLT1 #/BLT2
  #YLT <<- cbind(scale(allmean),scale(scores))%*%BLT
  YLT <<- cbind(as.matrix(allmean),as.matrix(scores))%*%BLT
  rnames <- rownames(YLT)
  thetaLT <- solve(MMVG)%*%(MMVP)%*%BLT
  thetaNLT <- thetaLT%*%(1/sqrt(t(thetaLT)%*%thetaLT))
  vGb <- t(thetaNLT)%*%MMVG%*%BLT
  #vGb <- t(theta)%*%MMVG%*%BLT
  VIDS <- t(BLT)%*%MMVP%*%BLT
  VDCG <- t(thetaNLT)%*%MMVG%*%thetaNLT
  #VDCG <- t(theta)%*%MMVG%*%theta
  vGv <- sqrt(abs(VDCG))
  bPb <- sqrt(abs(VIDS))
  corrAnalysis <- min(0.9999,vGb/(vGv%*%bPb))  
  selected <- YLT>=quantile(YLT,1-selval/100)
  selYLT <- YLT[selected]
  ord <- order(selYLT,decreasing=TRUE)
  selYLT <- selYLT[ord]
  selentry <-  allmean[selected,]
  nom.selentry <- rownames(allmean)[selected]
  selentry <- data.frame(as.matrix(selentry)[ord,])
  dimnames(selentry) <- list(nom.selentry[ord],nom.traits)
  MLTSI <- apply(selentry,2,mean)
  MLTall <- apply(allmean,2,mean)
  GAINbyMEANS=MLTSI-MLTall
  ks <-(100/selval)*(1/sqrt(2*pi))*exp((-qnorm(1-(selval/100))^2)/2)
  #ks <- if(selval==5) 2.063 else if(selval==10) 1.755 else 1.400
  zero<-rep(0, nrow(MVG))
   oaov <- c(theta1,zero)
   beta<-solve(MMVP)%*%MMVG%*%oaov
  Rsel<-PlotResp(t(beta)%*%MMVP%*%beta,"LT",selval)
  #Rsel<-ks*sqrt(t(beta)%*%MMVP%*%beta)
  H2<-vGb/(vGv*bPb)
  LTGain <- as.vector(ks*(MMVG%*%BLT)/as.numeric(sqrt(t(BLT)%*%MMVP%*%BLT)))
  selentry2 <- rbind(selentry,MLTSI,MLTall,GAINbyMEANS,LTGain)
  selLT <- data.frame(selentry2,c(selYLT,rep(NA,4)))
  colnames(selLT)[ncol(selLT)] <- "LT index"
  rnames <- c("Mean of Selected Individuals", "Mean of all Individuals",
            'Selection Differential',paste("Expected Genetic Gain for",paste(selval,"%",sep=""),sep=" "))
  rownames(selLT)[(length(selYLT)+1):(length(selYLT)+4)] <- rnames
  allLT <- data.frame(allmean,YLT)
   colnames(allLT)[ncol(allLT)] <- "LT index"
  cat("LANDE & THOMPSON SELECTION INDEX METHOD","\n",file=out)
  cat("\n",paste("GENETIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(MVG,2),file=out,col.names=T,append=T)
  cat("\n",paste("PHENOTYPIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(MVP,2),file=out,col.names=T,append=T)
  cat("\n",paste("MOLECULAR",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(MVM,2),file=out,col.names=T,append=T)
  if (corr != TRUE) {
  cat("\n\n","COVARIANCE BETWEEN THE LT SELECTION INDEX AND THE BREEDING VALUE: ", 
      round(vGb,3),"\n",file=out,append=T)
  cat("\n","VARIANCE OF THE LT SELECTION INDEX:                                 ", 
      round(VIDS,3),"\n",file=out,append=T)
  cat("\n","VARIANCE OF THE BREEDING VALUE:                                     ", 
      round(VDCG,3),"\n",file=out,append=T)
  cat("\n","CORRELATION BETWEEN THE LT SELECTION INDEX AND THE BREEDING VALUE:  ", 
      round(corrAnalysis,3),"\n",file=out,append=T)
  cat("\n", "RESPONSE TO SELECTION:                                             ", 
      round(Rsel,3), "\n", file = out, append = T)
  cat("\n", "HERITABILITY:                                                      ", 
      round(H2,3), "\n", file = out, append = T)
  }	  
  cat("\n\n",
      paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE LT SELECTION, MEANS AND GAINS FOR", paste(selval,"%",
            sep =""),sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(selLT,2),file=out,col.names=T,append=T)
  cat("\n\n","VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE LT SELECTION INDEX","\n",file=out,append=T)
  print.char.matrix(round(allLT,2),file=out,col.names=T,append=T)
  #file.show(out)
  write.csv(selentry,na="",file=outcsv)
  #write.csv(allLT,file=paste("all",outcsv,sep=""))
  allprint=allLT
  if(rawdata == TRUE){SummaryStats(traits, design, softR, one.env, block.ex, covG, covP,allprint, "BLUE_BLUP_H2", nom.traits)}
  save(MVG,MVP,vGb,VIDS,VDCG,corrAnalysis,Rsel,H2,MVM,selentry,allprint,file="forReport.RData")
  #cat("\n",paste("Output saved in the files:",paste(out,",",sep=""),outcsv,"and",paste("all",outcsv,sep=""),sep=" "),"\n")
  #rm(fcol, ntraits, mat, pos = ".GlobalEnv")  
}

