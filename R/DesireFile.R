DesireFile<-function(file.dat = NULL, file.wgt = NULL, selval = 5, design = "lattice", corr = FALSE, rawdata = TRUE,one.env=TRUE, block.ex, softR, file.covG, dir_root){

    traits <<- Read(file.dat, file.wgt,design,mis =  c(".","NA","na","N/A"), one.env, genomic=FALSE, markers=NULL, block.ex, softR)
	nom.traits <- colnames(traits)[-(1:fcol - 1)]
    
	n<-length(theta)
        covG <<- gCovariances(traits, design, softR,one.env,block.ex) 
        dimnames(covG) <- list(nom.traits, nom.traits)
        allmean <<- Means(traits,design, softR,one.env,block.ex)
        colnames(allmean) <- nom.traits
        covP <<- fCovariances(allmean)
        dimnames(covP) <- dimnames(covG)

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
########################
########################
#  write desire file
rGP <- covG
rGP[upper.tri(rGP)] <- t(covP)[upper.tri(covP)] # fill the upper triangular

covsT=list()
predictions1=list()
iter=0
for (k in fcol:(fcol+nrow(covG)-1)){
  iter=iter+1
  covsT[[iter]]=Adjcovs(k, data, design, softR, one.env, block.ex)
  predictions1[[iter]]=AdjBlup(k, data, design, softR, one.env, block.ex)
}


h2s=list()
predictions1=as.data.frame(do.call(cbind, predictions1))
for (p in 1:nrow(covG)){
  #h2S<-vc.g/(vc.g + (vc.e))
  h2s[[p]]<-covsT[[p]][1]/(covsT[[p]][1]+covsT[[p]][2])  
}
ew<- sqrt(apply(predictions1,2,var,na.rm=TRUE))/apply(predictions1,2,mean,na.rm=T)
#ew[-c(1:3)] <- ew[-c(1:3)] * -1

s1 <- data.frame("  Input file for program Desire."); colnames(s1)<-"m"
s2 <- data.frame("  Number of traits ...  "); colnames(s2)<-"m"
s3 <- data.frame(as.character(nrow(covG))); colnames(s3) <- "m"
s4 <- data.frame("  Trait names ..."); colnames(s4)<-"m"
s5 <- data.frame(as.character(colnames(covG))); colnames(s5)<-"m"
s6 <- data.frame("  Starting economic weights ..."); colnames(s6)<-"m"
s7 <- data.frame(as.character(ew)); colnames(s7)<-"m" 
s8 <- data.frame("  Heritabilities ..."); colnames(s8)<-"m"
s9 <- data.frame(as.character(round(as.numeric(unlist(h2s)),3))); colnames(s9)<-"m"
s10 <- data.frame("  Standard Deviations ..."); colnames(s10)<-"m"
s11 <- data.frame(as.character(rep(1,nrow(covG)))); colnames(s11)<-"m"
s12 <- data.frame("  Correlation matrix (rp on upper diagonal, rg on lower)."); colnames(s12)<-"m"
s13 <- apply(round(rGP,3),2,as.character)

s.1.12 <- rbind(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
s.1.12plus <- cbind(s.1.12,matrix("",ncol=ncol(s13) - 1, nrow=nrow(s.1.12)))
colnames(s.1.12plus) <- colnames(s13)
desireResult <- rbind(s.1.12plus,s13)

write.table(desireResult, file =file.path(dir_root,"Programs","Rcode","desire",paste0("desireResultOutSD.txt")),
            sep = "\t",row.names = FALSE, quote = FALSE, col.names = FALSE)
}