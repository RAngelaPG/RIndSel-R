PPGLPSI <-function (file.dat = NULL, file.wgt = NULL, selval = 10, design = "lattice", 
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
    
	cat("\n", "   Wait ...", "\n\n")    
	theta <- wgts[match(nom.traits,wgts[,1]),2]
	indrest <- wgts[match(nom.traits,wgts[,1]),3]	
	names(theta)<-nom.traits
	n<-length(theta)	    
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
    covGV <- MVG %*% theta
    Bs <- IcovP %*% covGV

    cone<-indrest[which(indrest!=0)]
    W<- diag(indrest)[, which(indrest > 0 | indrest < 0)]
    W[which(W!=0)]<-1
    cr<-diag(cone[length(cone)], length(cone)-1)
    minusc<- -cone[1:length(cone)-1]
    Cprim<-rbind(diag(cone[length(cone)], length(cone)-1),-cone[1:length(cone)-1])
    C <- t(MVG) %*% W
    Port <- IcovP %*% C %*% Cprim %*% solve(t(Cprim) %*% t(C) %*% IcovP %*% C %*% Cprim) %*% t(Cprim) %*% t(C)
    NR <- nrow(Port)
    A2 <- diag(NR)
    A <- A2 - Port
    #print(Port)
    #print(A)
    BPPGnn <- A %*% Bs
    BPPG <- BPPGnn
    YPPG <<- as.matrix(allmean) %*% BPPG
    rnames <- rownames(YPPG)
    vGb <- t(theta) %*% MVG %*% BPPG
    VIDS <- t(BPPG) %*% MVP %*% BPPG
    VDCG <- t(theta) %*% MVG %*% theta
    vGv <- sqrt(abs(VDCG))
    bPb <- sqrt(abs(VIDS))
    corrAnalysis <- min(0.9999, vGb/(vGv %*% bPb))
    H2<- (t(BPPG) %*% MVG %*% BPPG)/(t(BPPG) %*% MVP %*% BPPG)
    selected <- YPPG >= quantile(YPPG, 1 - selval/100)
    selYPPG <- YPPG[selected]
    ord <- order(selYPPG, decreasing = TRUE)
    selYPPG <- selYPPG[ord]
    selentry <- allmean[selected, ]
    nom.selentry <- rownames(allmean)[selected]
    selentry <- data.frame(as.matrix(selentry)[ord, ])
    dimnames(selentry) <- list(nom.selentry[ord], nom.traits)
    MPPGSI <- apply(selentry, 2, mean)
    MPPGall <- apply(allmean, 2, mean)
    GAINbyMEANS = MPPGSI - MPPGall
    ks <-(100/selval)*(1/sqrt(2*pi))*exp((-qnorm(1-(selval/100))^2)/2)
    PPGGain <- as.vector(ks * (covG %*% BPPG)/as.numeric(sqrt(t(BPPG) %*% covP %*% BPPG)))
	Rsel<-PlotResp(VIDS,"PPG",selval)
    #RS <- ks*sqrt(VIDS)		
    selentry2 <- rbind(selentry, MPPGSI, MPPGall, GAINbyMEANS,PPGGain)
    selPPG <- data.frame(selentry2, c(selYPPG, rep(NA, 4)))
    colnames(selPPG)[ncol(selPPG)] <- "PPG index"
    rnames <- c("Mean of Selected Individuals", "Mean of all Individuals", 
        "Selection Differential", paste("Expected Genetic Gain for", 
            paste(selval, "%", sep = ""), sep = " "))
    rownames(selPPG)[(length(selYPPG) + 1):(length(selYPPG) + 4)] <- rnames
    allPPG <- data.frame(allmean, YPPG)
    colnames(allPPG)[ncol(allPPG)] <- "PPG index"
    cat("THE KEMPTHORNE AND NORDSKOG (PPG) RESTRICTIVE SELECTION INDEX", 
        "\n", file = out)
    cat("\n", paste("GENETIC", mat.name, "MATRIX", sep = " "), 
        "\n", file = out, append = T)
    print.char.matrix(round(MVG, 2), file = out, col.names = T, 
        append = T)
    cat("\n", paste("PHENOTYPIC", mat.name, "MATRIX", sep = " "), 
        "\n", file = out, append = T)
    print.char.matrix(round(MVP, 2), file = out, col.names = T, 
        append = T)
    cat("\n\n", "COVARIANCE BETWEEN THE PPG SELECTION INDEX AND THE BREEDING VALUE: ", 
        round(vGb,3), "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE PPG SELECTION INDEX:                                 ", 
        round(VIDS,3), "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE BREEDING VALUE:                                      ", 
        round(VDCG,3), "\n", file = out, append = T)
    cat("\n", "CORRELATION BETWEEN THE PPG SELECTION INDEX AND THE BREEDING VALUE:  ", 
        round(corrAnalysis,3), "\n", file = out, append = T)
    cat("\n", "RESPONSE TO SELECTION:                                               ", 
        round(Rsel,3), "\n", file = out, append = T)
    cat("\n", "HERITABILITY:                                                        ", 
        round(H2,3), "\n", file = out, append = T)
    cat("\n\n", paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE PPG SELECTION, MEANS AND GAINS FOR", 
        paste(selval, "%", sep = ""), sep = " "), "\n", file = out, 
        append = T)
    print.char.matrix(round(selPPG, 2), file = out, col.names = T, 
        append = T)
    cat("\n\n", "VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE PPG SELECTION INDEX", 
        "\n", file = out, append = T)
    print.char.matrix(round(allPPG, 2), file = out, col.names = T, 
        append = T)
    #file.show(out)
    write.csv(selentry, na = "", file = outcsv)
    #write.csv(allPPG, file = paste("all", outcsv, sep = ""))
	allprint=allPPG
	if(rawdata == TRUE){SummaryStats(traits, design, softR, one.env, block.ex, covG, covP,allprint, "BLUE_BLUP_H2", nom.traits)}
	save(MVG,MVP,vGb,VIDS,VDCG,corrAnalysis,Rsel,H2,selentry,allprint,file="forReport.RData")
#    cat("\n", paste("Output saved in the files:", paste(out, 
#        ",", sep = ""), outcsv, "and", paste("all", outcsv, sep = ""), 
#        sep = " "), "\n")
    #rm(fcol, ntraits, mat, pos = ".GlobalEnv")	
}
