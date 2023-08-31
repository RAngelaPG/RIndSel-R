AdjBlup <-function (i, data, design=design, softR=softR, one.env=one.env, block.ex=block.ex){
####################################################INDIVIDUAL##########################################################################################################################################
####ACBD
if(one.env==T){  
  if(design=="ACBD"){
	if(softR=="AsremlR"){
		if(block.ex==T){fit<-asreml(as.formula(paste(names(data)[i]," ~ Block")), random=~ ENTRY + R + C, residual=~ ar1(R):ar1(C), data=data)}
		if(block.ex==F){fit<-asreml(as.formula(paste(names(data)[i]," ~ 1")), random=~ ENTRY + R + C, residual=~ ar1(R):ar1(C), data=data)}
	}else{
		if(block.ex==T){fit<-remlf90(as.formula(paste(names(data)[i]," ~ Block")), random=~ ENTRY + R + C, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
		if(block.ex==F){fit<-remlf90(as.formula(paste(names(data)[i]," ~ 1")), random=~ ENTRY + R + C, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
		}
	}
####Lattice
  if(design=="lattice"){
        fit <- lmer(data[, i] ~ 1 + (1 | REP/Block) + (1 | ENTRY), data = data, REML = TRUE)
    }
####RCBD
  if (design=="rcb"){	
        fit <- lmer(data[, i] ~ 1 + (1 | REP) + (1 | ENTRY), data = data, REML = TRUE)
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
	}else{
		if(block.ex==T){fit<-remlf90(as.formula(paste(names(data)[i]," ~ 1")), random=~ ENTRY + LOC + LxE + LxR + LxC + LxB, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
		if(block.ex==F){fit<-remlf90(as.formula(paste(names(data)[i]," ~ 1")), random=~ ENTRY + LOC + LxE + LxR + LxC, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
	}
  }
####Lattice
  if(fcol==5 && design=="lattice"){
      fit <- lmer(data[, i] ~ 1 + (1 | LOC/REP/Block) + (1 | LOC:ENTRY) + (1 | ENTRY), data = data, REML = TRUE)
    }
####RCBD
  if (fcol==4 && design=="rcb"){
	  fit <-lmer(data[, i] ~ 1 + (1 | LOC/REP) + (1 | LOC:ENTRY)  + (1 | ENTRY), data = data, REML = TRUE)
    }
####################################################COMBINED##########################################################################################################################################
}

if (design=="ACBD"){
		if(softR=="AsremlR"){
			means <- predict(fit,"ENTRY")$pvals[,2]
		}else{
			if(block.ex==F){
				mu=as.data.frame(fit$fixed$Intercept)$value
				means<-mu+as.data.frame(fit$ranef$ENTRY)$value
			}
			if(block.ex==T){
				mu=mean(as.data.frame(fit$fixed)[,1])
				means<-mu+as.data.frame(fit$ranef$ENTRY)$value
			}
		}
	}else{
		means<-fixef(fit)+ranef(fit)$ENTRY
	}  

    return(means)
}
