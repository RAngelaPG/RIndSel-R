Adjcovs <-function (i, data, design=design, softR=softR, one.env=one.env, block.ex=block.ex){
####################################################INDIVIDUAL##########################################################################################################################################
####ACBD
if(one.env==T){  
  if(design=="ACBD"){
	if(softR=="AsremlR"){
		if(block.ex==T){fit<-asreml(as.formula(paste(names(data)[i]," ~ Block")), random=~ ENTRY + R + C, residual=~ ar1(R):ar1(C), data=data)}
		if(block.ex==F){fit<-asreml(as.formula(paste(names(data)[i]," ~ 1")), random=~ ENTRY + R + C, residual=~ ar1(R):ar1(C), data=data)}
		cov.names <- c("Entry", "Residual","R","C")
		covs<-c(summary(fit)$var["ENTRY",1],summary(fit)$var["R:C!R",1],summary(fit)$var["R",1],summary(fit)$var["C",1])
	}else{
		if(block.ex==T){fit<-remlf90(as.formula(paste(names(data)[i]," ~ Block")), random=~ ENTRY + R + C, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
		if(block.ex==F){fit<-remlf90(as.formula(paste(names(data)[i]," ~ 1")), random=~ ENTRY + R + C, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
		cov.names <- c("Entry", "Residual","C")
		covs<-c(summary(fit)$var["ENTRY",1],summary(fit)$var["Residual",1],summary(fit)$var["C",1])
		}
	}
####Lattice
  if(design=="lattice"){
        varcorr <- VarCorr(lmer(data[, i] ~ 1 + (1 | REP/Block) + (1 | ENTRY), data = data, REML = TRUE))
        cov.names <- c("Entry", "Residual","Block (Rep)", "Rep")
        covs <- c(varcorr$ENTRY, attr(varcorr, "sc")^2,varcorr$'Block:REP', varcorr$REP)
    }
####RCBD
  if (design=="rcb"){	
        varcorr <- VarCorr(lmer(data[, i] ~ 1 + (1 | REP) + (1 | ENTRY), data = data, REML = TRUE))
        cov.names <- c("Entry", "Residual","Rep")
        covs <- c(varcorr$ENTRY, attr(varcorr, "sc")^2,varcorr$REP)
    }
####################################################INDIVIDUAL##########################################################################################################################################
}
####################################################COMBINED##########################################################################################################################################
####ACBD
if(one.env==F){
  if(design=="ACBD"){
	if(softR=="AsremlR"){
		if(block.ex==T){fit<-asreml(as.formula(paste(names(data)[i]," ~ 1")), random=~ ENTRY + LOC + LOC:ENTRY + LOC:R + LOC:C + LOC:Block, residual=~ LOC:ar1(R):ar1(C), data=data)}
		if(block.ex==F){fit<-asreml(as.formula(paste(names(data)[i]," ~ 1")), random=~ ENTRY + LOC + LOC:ENTRY + LOC:R + LOC:C, residual=~ LOC:ar1(R):ar1(C), data=data)}
		cov.names <- c("Entry", "Residual","LOC:ENTRY")
		covs<-c(summary(fit)$var["ENTRY",1],summary(fit)$var["LOC:R:C!R",1],summary(fit)$var["LOC:ENTRY",1],summary(fit)$var["LOC:C",1])
	}else{
		if(block.ex==T){fit<-remlf90(as.formula(paste(names(data)[i]," ~ 1")), random=~ ENTRY + LOC + LxE + LxR + LxC + LxB, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
		if(block.ex==F){fit<-remlf90(as.formula(paste(names(data)[i]," ~ 1")), random=~ ENTRY + LOC + LxE + LxR + LxC, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
		cov.names <- c("Entry","Residual","LOC:ENTRY")
		covs<-c(summary(fit)$var["ENTRY",1],summary(fit)$var["Residual",1],summary(fit)$var["LxE",1])
		}
  }
####Lattice
  if(fcol==5 && design=="lattice"){
      varcorr <- VarCorr(lmer(data[, i] ~ 1 + (1 | LOC/REP/Block) + (1 | LOC:ENTRY) + (1 | ENTRY), data = data, REML = TRUE))
      cov.names <- c("Entry", "Residual", "LOC:ENTRY")
      covs <- c(varcorr$ENTRY, attr(varcorr, "sc")^2, varcorr$'LOC:ENTRY')
    }
####RCBD
  if (fcol==4 && design=="rcb"){
	  varcorr <- VarCorr(lmer(data[, i] ~ 1 + (1 | LOC/REP) + (1 | LOC:ENTRY)  + (1 | ENTRY), data = data, REML = TRUE))
      cov.names <- c("Entry","Residual","LOC:ENTRY")
      covs <- c(varcorr$ENTRY, attr(varcorr, "sc")^2, varcorr$'LOC:ENTRY')
    }
####################################################COMBINED##########################################################################################################################################
}
    names(covs) <- cov.names
    return(covs)
}



