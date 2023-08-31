TallisIndex <-function (file.dat = NULL, file.wgt = NULL, selval = 5, design = "lattice", 
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
    theta <- Readwgt()
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
    BTallis1 <- as.vector(solve(MVP) %*% MVG %*% theta)
    BTallis2 <- as.vector(sqrt(t(BTallis1) %*% BTallis1))
    BTallis <- BTallis1/BTallis2
    YTallis <<- scale(allmean) %*% BTallis
    rnames <- rownames(YTallis)
    vGb <- t(theta) %*% MVG %*% BTallis
    VIDS <- t(BTallis) %*% MVP %*% BTallis
    VDCG <- t(theta) %*% MVG %*% theta
    vGv <- sqrt(abs(VDCG))
    bPb <- sqrt(abs(VIDS))
    corrAnalysis <- min(0.9999, vGb/(vGv %*% bPb))
    selected <- YTallis >= quantile(YTallis, 1 - selval/100)
    selYTallis <- YTallis[selected]
    ord <- order(selYTallis, decreasing = TRUE)
    selYTallis <- selYTallis[ord]
    selentry <- allmean[selected, ]
    nom.selentry <- rownames(allmean)[selected]
    selentry <- data.frame(as.matrix(selentry)[ord, ])
    dimnames(selentry) <- list(nom.selentry[ord], nom.traits)
    MTallisSI <- apply(selentry, 2, mean)
    MTallisall <- apply(allmean, 2, mean)
    GAINbyMEANS = MTallisSI - MTallisall
    ks <-(100/selval)*(1/sqrt(2*pi))*exp((-qnorm(1-(selval/100))^2)/2)
  #ks <- if(selval==5) 2.063 else if(selval==10) 1.755 else 1.400
    TallisGain = as.vector(ks * (covG %*% BTallis)/as.numeric(sqrt(t(BTallis) %*% 
        covP %*% BTallis)))
    selentry2 <- rbind(selentry, MTallisSI, MTallisall, GAINbyMEANS, 
        TallisGain)
    selTallis <- data.frame(selentry2, c(selYTallis, rep(NA, 4)))
    colnames(selTallis)[ncol(selTallis)] <- "Tallis index"
    rnames <- c("Mean of Selected Individuals", "Mean of all Individuals", 
        "Selection Differential", paste("Expected Genetic Gain for", 
            paste(selval, "%", sep = ""), sep = " "))
    rownames(selTallis)[(length(selYTallis) + 1):(length(selYTallis) + 
        4)] <- rnames
    allTallis <- data.frame(allmean, YTallis)
    colnames(allTallis)[ncol(selTallis)] <- "Tallis index"
    cat("Tallis SELECTION INDEX METHOD", "\n", file = out)
    cat("\n", paste("GENETIC", mat.name, "MATRIX", sep = " "), 
        "\n", file = out, append = T)
    print.char.matrix(round(MVG, 2), file = out, col.names = T, 
        append = T)
    cat("\n", paste("PHENOTYPIC", mat.name, "MATRIX", sep = " "), 
        "\n", file = out, append = T)
    print.char.matrix(round(MVP, 2), file = out, col.names = T, 
        append = T)
    cat("\n\n", "COVARIANCE BETWEEN THE Tallis SELECTION INDEX AND THE BREEDING VALUE: ", 
        vGb, "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE Tallis SELECTION INDEX:                               ", 
        VIDS, "\n", file = out, append = T)
    cat("\n", "VARIANCE OF THE BREEDING VALUE:                                      ", 
        VDCG, "\n", file = out, append = T)
    cat("\n", "CORRELATION BETWEEN THE Tallis SELECTION INDEX AND THE BREEDING VALUE:", 
        corrAnalysis, "\n", file = out, append = T)
    cat("\n\n", paste("VALUES OF THE TRAITS FOR SELECTED INDIVIDUALS AND THE VALUE OF THE Tallis SELECTION, MEANS AND GAINS FOR", 
        paste(selval, "%", sep = ""), sep = " "), "\n", file = out, 
        append = T)
    print.char.matrix(round(selTallis, 2), file = out, col.names = T, 
        append = T)
    cat("\n\n", "VALUES OF THE TRAITS FOR ALL INDIVIDUALS AND THE VALUE OF THE Tallis SELECTION INDEX", 
        "\n", file = out, append = T)
    print.char.matrix(round(allTallis, 2), file = out, col.names = T, 
        append = T)
    #file.show(out)
    write.csv(selentry, na = "", file = outcsv)
    #write.csv(allTallis, file = paste("all", outcsv, sep = ""))
	allprint=allTallis
	if(rawdata == TRUE){SummaryStats(traits, design, softR, one.env, block.ex, covG, covP,allprint, "BLUE_BLUP_H2", nom.traits)}
	save(MVG,MVP,vGb,VIDS,VDCG,corrAnalysis,Rsel,H2,selentry,allprint,file="forReport.RData")
    #cat("\n", paste("Output saved in the files:", paste(out,",", sep = ""), outcsv, "and", paste("all", outcsv, sep = ""), sep = " "), "\n")
    #rm(ntraits, mat, pos = ".GlobalEnv")
}

