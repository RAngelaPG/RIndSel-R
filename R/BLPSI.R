BLPSI <-function (file.dat = NULL, file.wgt = NULL, selval = 5, design = "lattice", 
    corr = FALSE, out = "out.txt", outcsv = "out.csv", rawdata = TRUE,one.env=TRUE, block.ex, softR,file.covG){
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
	#theta <- Readwgt()
    cat("\n", "   Wait ...", "\n\n")    
	theta <- wgts[match(nom.traits,wgts[,1]),2]
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
    }else {
        mat.name <- "COVARIANCE"
        MVG <- covG
        MVP <- covP
    }	
	BSmith <-  theta 
	YSmith <<- as.matrix(allmean) %*% BSmith
    rnames <- rownames(YSmith)
    vGb <- t(theta) %*% MVG %*% BSmith
    VIDS <- t(theta) %*% MVG %*% BSmith
    VDCG <- t(theta) %*% MVG %*% theta
    vGv <- sqrt(abs(VDCG))
    bPb <- sqrt(abs(VIDS))
	H2=sqrt((t(BSmith) %*% MVG %*% BSmith)/(t(BSmith)%*% MVP %*% BSmith))
    corrAnalysis=sqrt(H2)    
    #corrAnalysis <- min(0.9999, vGb/(vGv %*% bPb))
    selected <- YSmith >= quantile(YSmith, 1 - selval/100,na.rm=T)
    selYSmith <- YSmith[selected]
    ord <- order(selYSmith, decreasing = TRUE)
    selYSmith <- selYSmith[ord]
    selentry <- allmean[selected, ]
    nom.selentry <- rownames(allmean)[selected]
    selentry <- data.frame(as.matrix(selentry)[ord, ])
    dimnames(selentry) <- list(nom.selentry[ord], nom.traits)
    MSmithSI <- apply(selentry, 2, mean,na.rm=T)
    MSmithall <- apply(allmean, 2, mean,na.rm=T)
    GAINbyMEANS = MSmithSI - MSmithall
    ks <-(100/selval)*(1/sqrt(2*pi))*exp((-qnorm(1-(selval/100))^2)/2)
	Rsel<-PlotResp(VIDS,"BLPSI",selval)	
	#Rsel<-ks*sqrt(VIDS)
	H2<-sqrt((t(BSmith)%*%MVG%*%BSmith)/(t(BSmith)%*%MVP%*%BSmith))
    SmithGain = as.vector(ks * (covG %*% BSmith)/as.numeric(sqrt(VIDS)))
    selentry2 <- rbind(selentry, MSmithSI, MSmithall, GAINbyMEANS,SmithGain)
    selSmith <- data.frame(selentry2, c(selYSmith, rep(NA, 4)))
    colnames(selSmith)[ncol(selSmith)] <- "BLPSI"
    rnames <- c("Mean of Selected Individuals", "Mean of all Individuals", 
        "Selection Differential", paste("Expected Genetic Gain for", 
            paste(selval, "%", sep = ""), sep = " "))
    rownames(selSmith)[(length(selYSmith) + 1):(length(selYSmith) + 4)] <- rnames
    allSmith <- data.frame(allmean, YSmith)
    colnames(allSmith)[ncol(selSmith)] <- "BPLSI"
    cat("BPLSI: Base Linear Phenotypic Selection Index", "\n", file = out)
    cat("\n", paste("GENETIC", mat.name, "MATRIX", sep = " "), "\n", file = out, append = T)
    print.char.matrix(round(MVG, 2), file = out, col.names = T,append = T)
    cat("\n", paste("PHENOTYPIC", mat.name, "MATRIX", sep = " "), "\n", file = out, append = T)
    print.char.matrix(round(MVP, 2), file = out, col.names = T,append = T)
	if (corr != TRUE) {
    cat("\n\n", "COVARIANCE BETWEEN THE BLPSI SELECTION INDEX AND THE BREEDING VALUE: ", 
        round(vGb,3), "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE BLPSI SELECTION INDEX:                                 ", 
        round(VIDS,3), "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE BREEDING VALUE:                                        ", 
        round(VDCG,3), "\n", file = out, append = T)
    cat("\n", "CORRELATION BETWEEN THE BLPSI SELECTION INDEX AND THE BREEDING VALUE:  ", 
        round(corrAnalysis,3), "\n", file = out, append = T)
    cat("\n", "RESPONSE TO SELECTION:                                                 ", 
        round(Rsel,3), "\n", file = out, append = T)
    cat("\n", "HERITABILITY:                                                          ", 
        round(H2,3), "\n", file = out, append = T)	
    }
    cat("\n\n", paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE BLPSI SELECTION, MEANS AND GAINS FOR", 
        paste(selval, "%", sep = ""), sep = " "), "\n", file = out,append = T)
    print.char.matrix(round(selSmith, 2), file = out, col.names = T,append = T)
    cat("\n\n", "VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE BLPSI SELECTION INDEX", 
        "\n", file = out, append = T)
    print.char.matrix(round(allSmith, 2), file = out, col.names = T, append = T)
    #file.show(out)
    write.csv(selentry, na = "", file = outcsv)
    #write.csv(allSmith, file = paste("all", outcsv, sep = ""))	
	allprint=allSmith
	if(rawdata == TRUE){SummaryStats(traits, design, softR, one.env, block.ex, covG, covP,allprint, "BLUE_BLUP_H2", nom.traits)}
	save(MVG,MVP,vGb,VIDS,VDCG,corrAnalysis,Rsel,H2,selentry,allprint,file="forReport.RData")
    #cat("\n", paste("Output saved in the files:", paste(out, ",", sep = ""), outcsv, "and", paste("all", outcsv, sep = ""),  sep = " "), "\n")
    #rm(fcol, ntraits, mat, pos = ".GlobalEnv")	
}

