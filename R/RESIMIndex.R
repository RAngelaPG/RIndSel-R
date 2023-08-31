RESIMIndex <-function (file.dat = NULL, file.wgt = NULL, selval = 5, design = "lattice", 
    corr = FALSE, out = "out.txt", outcsv = "out.csv", rawdata = TRUE,one.env=T, block.ex, softR,file.covG){ 
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
	   
    theta <- wgts[match(nom.traits,wgts[,1]),2]
	indrest <- wgts[match(nom.traits,wgts[,1]),3]	
	names(theta)<-nom.traits
	n<-length(theta)
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
    C <- t(MVG) %*% W
    Port <- IcovP %*% C %*% solve(t(C) %*% IcovP %*% C) %*% t(C)
    NR <- nrow(Port)
    QI <- diag(NR)
    A2 <- QI - Port
    A <- A2 %*% IcovP %*% MVG
    BRESIM <- abs(svd(A)$u[, 1]) * theta
	WWE <- t(svd(A)$u[, 1]) %*% MVP %*% solve(MVG) %*% (solve(MVG) %*% svd(A)$u[, 1])
    thetaRESIM <- as.vector((solve(MVG) %*% MVP %*% svd(A)$u[, 1])) / WWE
	YRESIM <<- as.matrix(allmean) %*% BRESIM
    rnames <- rownames(YRESIM)
	vGb <- t(thetaRESIM) %*% MVP %*% svd(A)$u[,1]
    VIDS <- t(BRESIM) %*% MVP %*% BRESIM
    VDCG <- t(thetaRESIM) %*% MVG %*% thetaRESIM
    vGv <- sqrt(abs(VDCG))
    bPb <- sqrt(abs(VIDS))
	RVsi=sqrt(VIDS) 
    H2=sqrt((t(BRESIM) %*% MVG %*% BRESIM)/(t(BRESIM) %*% MVP %*% BRESIM))
    corrAnalysis=sqrt(H2) 
    selected <- YRESIM >= quantile(YRESIM, 1 - selval/100)
    selYRESIM <- YRESIM[selected]
    ord <- order(selYRESIM, decreasing = TRUE)
    selYRESIM <- selYRESIM[ord]
    selentry <- allmean[selected, ]
    nom.selentry <- rownames(allmean)[selected]
    selentry <- data.frame(as.matrix(selentry)[ord, ])
    dimnames(selentry) <- list(nom.selentry[ord], nom.traits)
    MRESIMSI <- apply(selentry, 2, mean)
    MRESIMall <- apply(allmean, 2, mean)
    GAINbyMEANS = MRESIMSI - MRESIMall
    ks <-(100/selval)*(1/sqrt(2*pi))*exp((-qnorm(1-(selval/100))^2)/2)
	Rsel<-PlotResp(VIDS,"RESIM",selval)
    #Rsel<-ks*RVsi		
    #H2<-(t(BRESIM)%*%MVG%*%BRESIM)/(t(BRESIM)%*%MVP%*%BRESIM)  
    RESIMGain <- as.vector(ks * (covG %*% BRESIM)/as.numeric(sqrt(t(BRESIM) %*% covP %*% BRESIM)))
    selentry2 <- rbind(selentry, MRESIMSI, MRESIMall, GAINbyMEANS,RESIMGain)
    selRESIM <- data.frame(selentry2, c(selYRESIM, rep(NA, 4)))
    colnames(selRESIM)[ncol(selRESIM)] <- "RESIM index"
    rnames <- c("Mean of Selected Individuals", "Mean of all Individuals", 
        "Selection Differential", paste("Expected Genetic Gain for", 
            paste(selval, "%", sep = ""), sep = " "))
    rownames(selRESIM)[(length(selYRESIM) + 1):(length(selYRESIM) + 4)] <- rnames
    allRESIM <- data.frame(allmean, YRESIM)
    colnames(allRESIM)[ncol(allRESIM)] <- "RESIM index"
    cat("THE RESIM SELECTION INDEX", "\n", file = out)
    cat("\n", paste("GENETIC", mat.name, "MATRIX", sep = " "), "\n", file = out, append = T)
    print.char.matrix(round(MVG, 2), file = out, col.names = T, append = T)
    cat("\n", paste("PHENOTYPIC", mat.name, "MATRIX", sep = " "), "\n", file = out, append = T)
    print.char.matrix(round(MVP, 2), file = out, col.names = T, append = T)
	if (corr != TRUE) {
    cat("\n\n", "COVARIANCE BETWEEN THE RESIM SELECTION INDEX AND THE BREEDING VALUE: ", 
        round(vGb,3), "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE RESIM SELECTION INDEX:                                 ", 
        round(VIDS,3), "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE BREEDING VALUE:                                        ", 
        round(VDCG,3), "\n", file = out, append = T)
    cat("\n", "CORRELATION BETWEEN THE RESIM SELECTION INDEX AND THE BREEDING VALUE:  ", 
        round(corrAnalysis,3), "\n", file = out, append = T)
	cat("\n", "RESPONSE TO SELECTION:                                                 ", 
	    round(Rsel,3), "\n", file = out, append = T)
    cat("\n", "HERITABILITY:                                                          ", 
	    round(H2,3), "\n", file = out, append = T)	
    }		
    cat("\n\n", paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE RESIM SELECTION, MEANS AND GAINS FOR", 
        paste(selval, "%", sep = ""), sep = " "), "\n", file = out, append = T)
    print.char.matrix(round(selRESIM, 2), file = out, col.names = T, append = T)
    cat("\n\n", "VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE RESIM SELECTION INDEX", 
        "\n", file = out, append = T)
    print.char.matrix(round(allRESIM, 2), file = out, col.names = T, append = T)
    #file.show(out)
    write.csv(selentry, na = "", file = outcsv)
    #write.csv(allRESIM, file = paste("all", outcsv, sep = ""))
	allprint=allRESIM
	if(rawdata == TRUE){SummaryStats(traits, design, softR, one.env, block.ex, covG, covP,allprint, "BLUE_BLUP_H2", nom.traits)}
	save(MVG,MVP,vGb,VIDS,VDCG,corrAnalysis,Rsel,H2,selentry,allprint,file="forReport.RData")
    #cat("\n", paste("Output saved in the files:", paste(out,",", sep = ""), outcsv, "and", paste("all", outcsv, sep = ""),sep = " "), "\n")
    #rm(fcol, ntraits, mat, pos = ".GlobalEnv")	
}

