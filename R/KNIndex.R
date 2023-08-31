KNIndex <-function (file.dat = NULL, file.wgt = NULL, selval = 5, design = "lattice", 
    corr = FALSE, out = "out.txt", outcsv = "out.csv", rawdata = TRUE,one.env=T, block.ex, softR, file.covG){
    if (rawdata == TRUE) {
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
	
    if (is.null(traits)) {
    traits <<- Read(file.dat, file.wgt,design,mis =  c(".","NA","na","N/A"), one.env, genomic=FALSE, markers=NULL, block.ex, softR)
	nom.traits <- colnames(traits)[-(1:fcol - 1)]    
    }
    #thetaKN <- Readwgt()
	#W <- Readrest(weights)    
	thetaKN <- wgts[match(nom.traits,wgts[,1]),2]
	indrest <- wgts[match(nom.traits,wgts[,1]),3]
    names(thetaKN)<-nom.traits	
	n<-length(thetaKN)
	W <- diag(n)[, which(indrest == 1)]
    cat("\n", "   Wait ...", "\n\n")
    if (is.null(covG) & is.null(covP) & is.null(allmean)) {
        covG <<- gCovariances(traits, design, softR,one.env,block.ex) 
        dimnames(covG) <- list(nom.traits, nom.traits)
        allmean <<- Means(traits,design, softR,one.env,block.ex)
        colnames(allmean) <- nom.traits
        covP <<- fCovariances(allmean)
        dimnames(covP) <- dimnames(covG)
    }

    if (corr == TRUE) {
        mat.name <- "CORRELATION"
	  MVGS <<- cov2cor(covG)
        MVG <- MVGS
        MVPS <<- cov2cor(covP)
        MVP <- MVPS   
 	}

    else {
        mat.name <- "COVARIANCE"
        MVG <- covG
        MVP <- covP
    }
    IcovP <- solve(MVP)
    covGV <- MVG %*% thetaKN
    Bs <- IcovP %*% covGV
    C <- t(MVG) %*% W
    Port <- IcovP %*% C %*% solve(t(C) %*% IcovP %*% C) %*% t(C)
    NR <- nrow(Port)
    A2 <- diag(NR)
    A <- A2 - Port
    BKN <- A %*% Bs
	YKN <<- as.matrix(allmean) %*% BKN
    rnames <- rownames(YKN)
    vGb <- t(thetaKN) %*% MVG %*% BKN
    VIDS <- t(thetaKN) %*% MVG %*% BKN
    VDCG <- t(thetaKN) %*% MVG %*% thetaKN
    vGv <- sqrt(abs(VDCG))
    bPb <- sqrt(abs(VIDS))
	RVsi=sqrt(VIDS)
    H2=sqrt((t(thetaKN) %*% MVG %*% thetaKN)/(t(thetaKN)%*% MVP %*% thetaKN))
    corrAnalysis=sqrt(H2)
    selected <- YKN >= quantile(YKN, 1 - selval/100)
    selYKN <- YKN[selected]
    ord <- order(selYKN, decreasing = TRUE)
    selYKN <- selYKN[ord]
    selentry <- allmean[selected, ]
    nom.selentry <- rownames(allmean)[selected]
    selentry <- data.frame(as.matrix(selentry)[ord, ])
    dimnames(selentry) <- list(nom.selentry[ord], nom.traits)
    MKNSI <- apply(selentry, 2, mean)
    MKNall <- apply(allmean, 2, mean)
    GAINbyMEANS = MKNSI - MKNall
    ks <-(100/selval)*(1/sqrt(2*pi))*exp((-qnorm(1-(selval/100))^2)/2)
	Rsel<-PlotResp(VIDS,"KN",selval)
    #Rsel<-ks*RVsi	
	KNGain <- as.vector(ks * (covG %*% BKN)/as.numeric(sqrt(t(BKN) %*% covP %*% BKN)))
    selentry2 <- rbind(selentry, MKNSI, MKNall, GAINbyMEANS, KNGain)
    selKN <- data.frame(selentry2, c(selYKN, rep(NA, 4)))
    colnames(selKN)[ncol(selKN)] <- "KN index"
    rnames <- c("Mean of Selected Individuals", "Mean of all Individuals", 
        "Selection Differential", paste("Expected Genetic Gain for", 
            paste(selval, "%", sep = ""), sep = " "))
    rownames(selKN)[(length(selYKN) + 1):(length(selYKN) + 4)] <- rnames
    allKN <- data.frame(allmean, YKN)
    colnames(allKN)[ncol(allKN)] <- "KN index"
    cat("THE KEMPTHORNE AND NORDSKOG (KN) RESTRICTIVE SELECTION INDEX", "\n", file = out)
    cat("\n", paste("GENETIC", mat.name, "MATRIX", sep = " "),"\n", file = out, append = T)
    print.char.matrix(round(MVG, 2), file = out, col.names = T,append = T)
    cat("\n", paste("PHENOTYPIC", mat.name, "MATRIX", sep = " "),"\n", file = out, append = T)
    print.char.matrix(round(MVP, 2), file = out, col.names = T, append = T)
	if (corr != TRUE) {
    cat("\n\n", "COVARIANCE BETWEEN THE KN SELECTION INDEX AND THE BREEDING VALUE: ", 
        round(vGb,3), "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE KN SELECTION INDEX:                                 ", 
        round(VIDS,3), "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE BREEDING VALUE:                                     ", 
        round(VDCG,3), "\n", file = out, append = T)
    cat("\n", "CORRELATION BETWEEN THE KN SELECTION INDEX AND THE BREEDING VALUE:  ", 
        round(corrAnalysis,3), "\n", file = out, append = T)
	cat("\n", "RESPONSE TO SELECTION:                                              ", 
        round(Rsel,3), "\n", file = out, append = T)
    cat("\n", "HERITABILITY:                                                       ", 
        round(H2,3), "\n", file = out, append = T)
    }		
    cat("\n\n", paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE KN SELECTION, MEANS AND GAINS FOR", 
        paste(selval, "%", sep = ""), sep = " "), "\n", file = out,append = T)
    print.char.matrix(round(selKN, 2), file = out, col.names = T,append = T)
    cat("\n\n", "VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE KN SELECTION INDEX","\n", file = out, append = T)
    print.char.matrix(round(allKN, 2), file = out, col.names = T, append = T)
    #file.show(out)
    write.csv(selentry, na = "", file = outcsv)
    #write.csv(allKN, file = paste("all", outcsv, sep = ""))
	allprint=allKN
	if(rawdata == TRUE){SummaryStats(traits, design, softR, one.env, block.ex, covG, covP,allprint, "BLUE_BLUP_H2", nom.traits)}
	save(MVG,MVP,vGb,VIDS,VDCG,corrAnalysis,Rsel,H2,selentry,allprint,file="forReport.RData")
    #cat("\n", paste("Output saved in the files:", paste(out,",", sep = ""), outcsv, "and", paste("all", outcsv, sep = ""),sep = " "), "\n")
    #rm(fcol, ntraits, mat, pos = ".GlobalEnv")		
}

