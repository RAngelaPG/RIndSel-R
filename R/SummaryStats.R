SummaryStats<-function(data, design, softR, one.env, block.ex, covG, covP,BLUEs, outcsv, nom.traits){

	BLUEsSc<-as.data.frame(cbind(levels(data$ENTRY),scale(BLUEs),Crossf=rep("X",dim(BLUEs)[1]),Type=rep("BLUE",dim(BLUEs)[1])))
	colnames(BLUEsSc)[1]="Entry"
	
	BLUEs<-as.data.frame(cbind(levels(data$ENTRY),BLUEs,Crossf=rep("X",dim(BLUEs)[1]),Type=rep("BLUE",dim(BLUEs)[1])))
	colnames(BLUEs)[1]="Entry"
	
	BLUPs <- MeansBP(data,design, softR,one.env,block.ex)
	colnames(BLUPs) <- nom.traits
	
	BLUPsSc<-as.data.frame(cbind(levels(data$ENTRY),scale(BLUPs),Crossf=rep("X",dim(BLUPs)[1]),Type=rep("BLUP",dim(BLUPs)[1])))
	colnames(BLUPsSc)[1]="Entry"
	
	BLUPs<-as.data.frame(cbind(levels(data$ENTRY),BLUPs,Crossf=rep("X",dim(BLUPs)[1]),Type=rep("BLUP",dim(BLUPs)[1])))
	colnames(BLUPs)[1]="Entry"
	
	covsT=list()
	iter=0
	for (k in fcol:(fcol+nrow(covG)-1)){
		iter=iter+1
		covsT[[iter]]=Adjcovs(k, data, design, softR, one.env, block.ex)
	}
	
	if (one.env==FALSE){
	nameloc=levels(data[,"LOC"])
	H2env=cbind(h2s=rep(0,(length(nameloc)*length(nom.traits))),trait=rep(nom.traits,length(nameloc)),EnvTrialf=rep(nameloc, each=length(nom.traits)))
	iter2=0
	for ( nl in 1:length(nameloc)){
		datatmp=subset(data, LOC==nameloc[nl])
		for (nt in fcol:(fcol+nrow(covG)-1)){
		iter2=iter2+1
		covEnv=Adjcovs(nt, datatmp, design, softR, TRUE, block.ex)
		#h2S<-vc.g/(vc.g + (vc.e))
		H2env[iter2,1]<-round(covEnv[1]/(covEnv[1]+covEnv[2]),4)
		}
	}
		h2s=list()
		for (p in 1:nrow(covG)){
			#h2S<-vc.g/(vc.g + (vc.ge/nenv)+(vc.e/nenv*nrep))
			if (design=="ACBD"){
				h2s[[p]]<-covsT[[p]][1]/(covsT[[p]][1]+ (covsT[[p]][3]/length(nameloc)) +(covsT[[p]][2]/length(nameloc)))
			}else{
				h2s[[p]]<-covsT[[p]][1]/(covsT[[p]][1]+ (covsT[[p]][3]/length(nameloc)) +(covsT[[p]][2]/(length(nameloc)*nlevels(data$REP))))
			}
		}
		h2s=as.data.frame(cbind(unlist(h2s),nom.traits))
		colnames(h2s)=c("H2","traits")
		l <- list("predictions" = BLUEs, "predictionsSc" = BLUEsSc, 
			"predictionsHyb"= BLUPs, "predictionsHybSc"=BLUPsSc,
			"H2"=h2s, "G"=covG, "P"=covP,"H2.envs"=H2env)
	
	}else{	
		h2s=list()
		for (p in 1:nrow(covG)){
			#h2S<-vc.g/(vc.g + (vc.e))
			h2s[[p]]<-covsT[[p]][1]/(covsT[[p]][1]+covsT[[p]][2])  
		}
		h2s=as.data.frame(cbind(unlist(h2s),nom.traits))
		colnames(h2s)=c("H2","traits")
		l <- list("predictions" = BLUEs, "predictionsSc" = BLUEsSc, 
			"predictionsHyb"= BLUPs, "predictionsHybSc"=BLUPsSc,
			"H2"=h2s, "G"=covG, "P"=covP)
	
	}
	
	write.xlsx(l, file=file.path(paste0("results",outcsv,".xlsx")))


}