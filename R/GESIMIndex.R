GESIMIndex <- function(file.dat=NULL,file.wgt=NULL,selval=5,design="lattice",corr=FALSE,method=c("vanraden","classic"),out="out.txt",outcsv="out.csv",rawdata=TRUE,markers=NULL,one.env=T, block.ex, softR, file.covG){
    if(rawdata==TRUE){
    traits <<- NULL
    traits2 <<- NULL
    covG <<- NULL
    covP <<- NULL
    allmean <<- NULL
    gebv <<- NULL
  }else{
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
	  markers_actualizados<<- ReadMarkers(markers) 
	  gebv <<- GEBV(P=covP,G=covG,X=markers_actualizados,method=method)
	  cov.ebv <<- gebv$ga 
	  dimnames(cov.ebv) <- dimnames(covG)
	  }
	
  if(is.null(traits)){
    traits <<- Read(file.dat, file.wgt,design,mis =  c(".","NA","na","N/A"), one.env, genomic=TRUE, markers=markers, block.ex, softR)
	nom.traits <- colnames(traits)[-(1:fcol - 1)]
  }
  
      theta <- Readsgn()
      cat("\n","   Wait ...","\n\n")
  if(is.null(covG)&is.null(covP)&is.null(allmean)&is.null(gebv)){
    covG <<- gCovariances(traits, design, softR,one.env,block.ex) 
    dimnames(covG) <- list(nom.traits, nom.traits)
    allmean <<- Means(traits,design, softR,one.env,block.ex)
    colnames(allmean) <- nom.traits
    covP <<- fCovariances(allmean)
    dimnames(covP) <- dimnames(covG)
    gebv <<- GEBV(P=covP,G=covG,X=markers_actualizados,method=method)
    cov.ebv <<- gebv$ga 
    dimnames(cov.ebv) <- dimnames(covG)
  }
  
  if(corr==TRUE) {
    mat.name <- "CORRELATION"
    MVGS <<- cov2cor(covG)
    MVG <- MVGS
    MVPS <<- cov2cor(covP)
    MVP <- MVPS
    MVMS <<- cov2cor(cov.ebv)
    MVM <- MVMS
    
  } else{
    mat.name <- "COVARIANCE"
    MVG <- covG
    MVP <- covP
    MVM <- cov.ebv     
  }
  
  si=as.vector(wgts[,3])
  sinna=si[!is.na(si)]
  EigVecGESIM <- IndCoef(theta=sinna,Gamma=MVM,G=MVG,evd=T)
  EigVecGESIM <- EigVecGESIM/as.vector(sqrt(t(EigVecGESIM)%*%EigVecGESIM))
  ebv <- gebv$ebv
  YGESIM <<- scale(ebv)%*%EigVecGESIM
  rnames <- rownames(YGESIM)
  thetaGESIM <- solve(MVM)%*%(MVP)%*%EigVecGESIM
  thetaNGESIM <- thetaGESIM%*%(1/sqrt(t(thetaGESIM)%*%thetaGESIM))
  vGb <- t(thetaNGESIM)%*%MVM%*%EigVecGESIM
  VIDS <- t(EigVecGESIM)%*%MVP%*%EigVecGESIM
  VDCG <- t(thetaNGESIM)%*%MVM%*%thetaNGESIM
  vGv <- sqrt(abs(VDCG))
  bPb <- sqrt(abs(VIDS))
  corrAnalysis <- min(0.9999,ifelse(corr!=TRUE,vGb/(vGv%*%bPb),vGb))  
  selected <- YGESIM>=quantile(YGESIM,1-selval/100)
  selYGESIM <- YGESIM[selected]
  ord <- order(selYGESIM,decreasing=TRUE)
  selYGESIM <- selYGESIM[ord]
  selentry <-  allmean[selected,]
  nom.selentry <- rownames(allmean)[selected]
  selentry <- data.frame(as.matrix(selentry)[ord,])
  dimnames(selentry) <- list(nom.selentry[ord],nom.traits)
  MGESIMSI <- apply(selentry,2,mean)
  MGESIMall <- apply(allmean,2,mean)
  GAINbyMEANS=MGESIMSI-MGESIMall  
  ks <-(100/selval)*(1/sqrt(2*pi))*exp((-qnorm(1-(selval/100))^2)/2)
  #ks <- if(selval==5) 2.063 else if(selval==10) 1.755 else 1.400
  GESIMGain <- as.vector(ks*(MVM%*%EigVecGESIM)/as.numeric(sqrt(VIDS)))
  selentry2 <- rbind(selentry,MGESIMSI,MGESIMall,GAINbyMEANS,GESIMGain)
  selGESIM <- data.frame(selentry2,c(selYGESIM,rep(NA,4)))
  colnames(selGESIM)[ncol(selGESIM)] <- "GESIM index"
  rnames <- c("Mean of Selected Individuals", "Mean of all Individuals",
            'Selection Differential',paste("Expected Genetic Gain for",paste(selval,"%",sep=""),sep=" "))
  rownames(selGESIM)[(length(selYGESIM)+1):(length(selYGESIM)+4)] <- rnames
  allGESIM <- data.frame(allmean,YGESIM)
   colnames(allGESIM)[ncol(allGESIM)] <- "GESIM index"
  cat("GESIM SELECTION INDEX METHOD","\n",file=out)
  cat("\n",paste("GENETIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(MVG,2),file=out,col.names=T,append=T)
  cat("\n",paste("PHENOTYPIC",mat.name,"MATRIX",sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(MVP,2),file=out,col.names=T,append=T)
  cat("\n","GEBV COVARIANCE MATRIX","\n",file=out,append=T)
  print.char.matrix(round(MVM,2),file=out,col.names=T,append=T)
  cat("\n\n","COVARIANCE BETWEEN THE GESIM SELECTION INDEX AND THE BREEDING VALUE: ", vGb,"\n",file=out,append=T)
  cat("\n","VARIANCE OF THE GESIM SELECTION INDEX:                                 ", VIDS,"\n",file=out,append=T)
  cat("\n","VARIANCE OF THE BREEDING VALUE:                                        ", VDCG,"\n",file=out,append=T)
  cat("\n","CORRELATION BETWEEN THE GESIM SELECTION INDEX AND THE BREEDING VALUE:  ", corrAnalysis,"\n",file=out,append=T)
  cat("\n\n",
      paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE GESIM SELECTION, MEANS AND GAINS FOR", paste(selval,"%",
            sep =""),sep=" "),"\n",file=out,append=T)
  print.char.matrix(round(selGESIM,2),file=out,col.names=T,append=T)
  cat("\n\n","VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE GESIM SELECTION INDEX","\n",file=out,append=T)
  print.char.matrix(round(allGESIM,2),file=out,col.names=T,append=T)
  #file.show(out)
  write.csv(selentry,na="",file=outcsv)
  #write.csv(allGESIM,file=paste("all",outcsv,sep=""))
  allprint=allGESIM
  if(rawdata == TRUE){SummaryStats(traits, design, softR, one.env, block.ex, covG, covP,allprint, "BLUE_BLUP_H2", nom.traits)}
  save(MVG,MVP,vGb,VIDS,VDCG,corrAnalysis,Rsel,H2,MVM,selentry,allprint,file="forReport.RData")
  #cat("\n",paste("Output saved in the files:",paste(out,",",sep=""),outcsv,"and",paste("all",outcsv,sep=""),sep=" "),"\n")
  #rm(mat, pos = ".GlobalEnv")
}

