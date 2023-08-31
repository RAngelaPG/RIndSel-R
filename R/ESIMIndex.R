ESIMIndex <-function (file.dat = NULL, file.wgt = NULL,markers=NULL, selval = 5, design = "lattice", 
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
	names(theta)<-nom.traits
	n<-length(theta)
    
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

    svdIPG <- svd(solve(MVP) %*% MVG)
    EigVecESIM <- as.matrix(abs(svdIPG$u[, 1])) * theta
	YESIM <<- as.matrix(allmean) %*% EigVecESIM
    rnames <- rownames(YESIM)
	WWE <- t(svdIPG$u[, 1]) %*% MVP %*% solve(MVG) %*% (solve(MVG) %*% svdIPG$u[, 1])
    thetaESIM <- as.vector((solve(MVG) %*% MVP %*% svdIPG$u[, 1])) / WWE
    vGb <- t(thetaESIM) %*% MVP %*% svdIPG$u[,1]
    VIDS <- t(EigVecESIM) %*% MVP %*% EigVecESIM
    VDCG <- t(thetaESIM) %*% MVG %*% thetaESIM
    vGv <- sqrt(abs(VDCG))
    bPb <- sqrt(abs(VIDS))	
	RVsi=sqrt(VIDS) 
    H2=sqrt((t(EigVecESIM) %*% MVG %*% EigVecESIM)/(t(EigVecESIM) %*% MVP %*% EigVecESIM))
    corrAnalysis=sqrt(H2)
    selected <- YESIM >= quantile(YESIM, 1 - selval/100)
    selYESIM <- YESIM[selected]
    ord <- order(selYESIM, decreasing = TRUE)
    selYESIM <- selYESIM[ord]
    selentry <- allmean[selected, ]	
    nom.selentry <- rownames(allmean)[selected]	
    selentry <- data.frame(as.matrix(selentry)[ord, ])	
    dimnames(selentry) <- list(nom.selentry[ord], nom.traits)
    MESIMSI <- apply(selentry, 2, mean)
    MESIMall <- apply(allmean, 2, mean)
    GAINbyMEANS = MESIMSI - MESIMall	
    ks <-(100/selval)*(1/sqrt(2*pi))*exp((-qnorm(1-(selval/100))^2)/2)
	Rsel<-PlotResp(VIDS,"ESIM",selval)
    #Rsel<-ks*sqrt(VIDS)
	ESIMGain <- as.vector(ks * (covG %*% EigVecESIM)/as.numeric(sqrt(t(EigVecESIM) %*% covP %*% EigVecESIM)))
    selentry2 <- rbind(selentry, MESIMSI, MESIMall, GAINbyMEANS,ESIMGain)
    selESIM <- data.frame(selentry2, c(selYESIM, rep(NA, 4)))
    colnames(selESIM)[ncol(selESIM)] <- "ESIM index"
    rnames <- c("Mean of Selected Individuals", "Mean of all Individuals", 
        "Selection Differential", paste("Expected Genetic Gain for", 
            paste(selval, "%", sep = ""), sep = " "))
    rownames(selESIM)[(length(selYESIM) + 1):(length(selYESIM) + 
        4)] <- rnames		
    allESIM <- data.frame(allmean, YESIM)
    colnames(allESIM)[ncol(allESIM)] <- "ESIM index"
    cat("ESIM SELECTION INDEX METHOD", "\n", file = out)
    cat("\n", paste("GENETIC", mat.name, "MATRIX", sep = " "), 
        "\n", file = out, append = T)
    print.char.matrix(round(MVG, 2), file = out, col.names = T, 
        append = T)
    cat("\n", paste("PHENOTYPIC", mat.name, "MATRIX", sep = " "), 
        "\n", file = out, append = T)
    print.char.matrix(round(MVP, 2), file = out, col.names = T, 
        append = T)
	if (corr != TRUE) {
    cat("\n\n", "COVARIANCE BETWEEN THE ESIM SELECTION INDEX AND THE BREEDING VALUE: ", 
        round(vGb,3), "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE ESIM SELECTION INDEX:                                 ", 
        round(VIDS,3), "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE BREEDING VALUE:                                       ", 
        round(VDCG,3), "\n", file = out, append = T)
    cat("\n", "CORRELATION BETWEEN THE ESIM SELECTION INDEX AND THE BREEDING VALUE:  ", 
        round(corrAnalysis,3), "\n", file = out, append = T)
	cat("\n", "RESPONSE TO SELECTION:                                                ", 
        round(Rsel,3), "\n", file = out, append = T)
    cat("\n", "HERITABILITY:                                                         ", 
        round(H2,3), "\n", file = out, append = T)
    }		
    cat("\n\n", paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE ESIM SELECTION, MEANS AND GAINS FOR", 
        paste(selval, "%", sep = ""), sep = " "), "\n", file = out, 
        append = T)
    print.char.matrix(round(selESIM, 2), file = out, col.names = T, 
        append = T)
    cat("\n\n", "VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE ESIM SELECTION INDEX", 
        "\n", file = out, append = T)
    print.char.matrix(round(allESIM, 2), file = out, col.names = T, 
        append = T)
    #file.show(out)
    write.csv(selentry, na = "", file = outcsv)
    #write.csv(allESIM, file = paste("all", outcsv, sep = ""))
	allprint=allESIM
	if(rawdata == TRUE){SummaryStats(traits, design, softR, one.env, block.ex, covG, covP,allprint, "BLUE_BLUP_H2", nom.traits)}
	save(MVG,MVP,vGb,VIDS,VDCG,corrAnalysis,Rsel,H2,selentry,allprint, file="forReport.RData")
    #cat("\n", paste("Output saved in the files:", paste(out,",", sep = ""), outcsv, "and", paste("all", outcsv, sep = ""),sep = " "), "\n")
    #rm(fcol, ntraits, mat, pos = ".GlobalEnv")		
}

