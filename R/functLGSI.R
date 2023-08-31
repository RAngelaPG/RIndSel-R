ReadTestingLGSI <-function(Tmarkers=NULL)
{
  fname=Tmarkers
  Tmarkers <- read.csv(file=fname,na.strings = c(".","NA"))
  T1<-Tmarkers[,1]
  recode2 <- function(i){
                cod2 <- ifelse(Tmarkers[,i]==0,-1,
                        ifelse(Tmarkers[,i]==1,0,
                        ifelse(Tmarkers[,i]==2,1,NA)))
    return(cod2)
  }
if(Tmarkers[3,1]>2) Tmarkers <- Tmarkers[,-1]
nmark <- ncol(Tmarkers)
if(any(Tmarkers==2, na.rm=TRUE)) {codes2 <- sapply(1:nmark,recode2)} 
   else {codes2 <- Tmarkers}
  colnames(codes2) <- paste("M",1:nmark,sep="")
if(any(is.na(codes2))==TRUE) {
         nas <- which(is.na(codes2),arr.ind=T)
         for(i in 1:nrow(nas)){
             codes2[nas[i,1],nas[i,2]] <- 0; 
         }
} else {codes2}  
  mcodes2 = as.matrix(cbind(T1, codes2))
  return(mcodes2)
}

GenRelLGSI <- function(markers,sel=c('f2','dh')){
  X <- markers
  XXT <- tcrossprod(X)
  m <- ncol(X)
  r <- nrow(X)
  p <- sqrt(sapply(1:ncol(X),function(i)table(X[,i])[1])/nrow(X)) 
  if(sel[1]=='f2'){
    Cma <- sum(2*p*(1-p))
    phi <- XXT*2/Cma      
   } else{
      Cma <- sum(4*p*(1-p))
      phi <- XXT*4/(3*Cma)-matrix(1,r,r)
    }
    return(list(Phi=phi, cma=Cma))
    }

GEBVLGSI <- function(P,G,X,sel="f2",method=c("vanraden","classic")){
  mph <- GenRelLGSI(X,sel=sel)
  tn <- nrow(G)
  g <- nrow(mph$Phi)
  if(method[1]=="vanraden"){
     cm <- mph$cma
  }else{
     cm<-1  
  }
  covgbv <- G%x%mph$Phi 
  R <- P-G
  L <- R%*%solve(G) 
  Vy <-covgbv+R%x%diag(g)  
  Z <- diag(tn)%x%X
  A <- Z%*%t(Z)
  B <- A%*%solve(A/cm+L%x%diag(g))/cm 
  y <- apply(allmean,2,scale, scale=F)
  T <- solve(Vy)-solve(Vy)%*%diag(tn*g)%*%solve(t(diag(tn*g))%*%solve(Vy)%*%diag(tn*g))%*%t(diag(tn*g))%*%solve(Vy)
  yq<-as.matrix(allmean, ncol=1)  
  ahat <- B%*%as.vector(as.matrix(y))
  ahat <- matrix(ahat,ncol=tn,byrow=FALSE)  
  uhat <- (1/cm*t(Z)%*%solve(A/cm+L%x%diag(g)))%*%as.vector(as.matrix(y))
  uhat <- matrix(uhat, ncol=tn, byrow=FALSE)
  vargamma<-1/g*(t(scale(ahat, scale=FALSE))%*%mph$Phi%*%scale(ahat, scale=FALSE))
  pev <- covgbv%*%solve(Vy)%*%covgbv
  covgebv <- ginv((diag(tn)%x%diag(g))+(G%*%solve(R))%x%mph$Phi)%*%(G%x%mph$Phi)
  congebv <- diag(tn*g)-ginv(covgbv)%*%pev 
  Gama <- matrix(0,tn,tn)
  for(i in 1:tn){
    for(j in i:tn){
      Gama[i,j] <- mean(diag(covgebv[(1:g)+(i-1)*g,(1:g)+(j-1)*g]))
    }    
  }  
  Gama <- t(Gama)+Gama
  diag(Gama) <- 0.5*diag(Gama)
  dimnames(Gama) <- dimnames(covG)
  dimnames(vargamma) <- dimnames(covG)  
  return(list(cov=covgebv,conf=congebv,pev=pev,ga=Gama,ebv=ahat, ebvs=uhat, Vga=vargamma, cm=cm))
}
