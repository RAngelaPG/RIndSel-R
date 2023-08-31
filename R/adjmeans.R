adjmeans <-function (i, data, design=design, softR=softR, one.env=one.env, block.ex=block.ex){ 
####################################################INDIVIDUAL##########################################################################################################################################
####ACBD
if(one.env==T){
  if(design=="ACBD"){
	if(softR=="AsremlR"){
		if(block.ex==T){fit<-asreml(as.formula(paste(names(data)[i]," ~ Block + ENTRY")), random=~ R + C, residual=~ ar1(R):ar1(C), data=data)}
		if(block.ex==F){fit<-asreml(as.formula(paste(names(data)[i]," ~ ENTRY")), random=~ R + C, residual=~ ar1(R):ar1(C), data=data)}
	}else{
		if(block.ex==T){fit<-remlf90(as.formula(paste(names(data)[i]," ~ Block + ENTRY")), random=~ R + C, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
		if(block.ex==F){fit<-remlf90(as.formula(paste(names(data)[i]," ~ ENTRY")), random=~ R + C, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
		}
	}
####Lattice	
  if(design=="lattice"){
     fit<-lmer(data[, i] ~ ENTRY + (1 | REP/Block), data = data, REML = TRUE)
	}
###RCBD
  if (design=="rcb"){
     fit<-lmer(data[, i] ~ ENTRY + (1 | REP), data = data, REML = TRUE)
	}
####################################################INDIVIDUAL##########################################################################################################################################
}

####################################################COMBINED##########################################################################################################################################
####ACBD
if(one.env==F){
  if(design=="ACBD"){
	if(softR=="AsremlR"){
		if(block.ex==T){fit<-asreml(as.formula(paste(names(data)[i]," ~ ENTRY")), random=~ LOC + LOC:ENTRY + LOC:R + LOC:C + LOC:Block, residual=~ LOC:ar1(R):ar1(C), data=data)}
		if(block.ex==F){fit<-asreml(as.formula(paste(names(data)[i]," ~ ENTRY")), random=~ LOC + LOC:ENTRY + LOC:R + LOC:C, residual=~ LOC:ar1(R):ar1(C), data=data)}
	}else{
		if(block.ex==T){fit<-remlf90(as.formula(paste(names(data)[i]," ~ ENTRY")), random=~ LOC + LxE + LxR + LxC + LxB, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
		if(block.ex==F){fit<-remlf90(as.formula(paste(names(data)[i]," ~ ENTRY")), random=~ LOC + LxE + LxR + LxC, spatial = list(model = 'AR',coord = data[, c('Row','Col')]), data= data, method="em")}
	 }
}
####Lattice	
  if(fcol==5 && design=="lattice"){
     fit<-lmer(data[, i] ~ ENTRY + (1 | LOC/REP/Block)+(1|LOC:ENTRY), data = data, REML = TRUE)
	 }
####RCBD
  if (fcol==4 && design=="rcb"){
      fit<-lmer(data[, i] ~ENTRY + (1 | LOC/REP)+(1|LOC:ENTRY), data = data, REML = TRUE)
	  }
####################################################COMBINED##########################################################################################################################################
}
	if (design=="ACBD"){
		if(softR=="AsremlR"){
			means <- predict(fit,"ENTRY")$pvals[,2]
		}else{
			means <- fixef(fit)$ENTRY
		}
	}else{
		means <- tapply(model.matrix(fit)%*%fixef(fit),data$ENTRY[!is.na(data[,i])],mean)
	}  
  return(means)
}

