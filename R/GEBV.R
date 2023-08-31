GEBV <- function(P,G,X,sel="f2",method=c("classic","vanraden")){
  Phi <- GenRel(markers=X,sel=sel)
  t <- nrow(G)
  g <- nrow(Phi)
  if(method[1]=="vanraden"){
   p <- sqrt(sapply(1:ncol(X),function(i)table(X[,i])[1])/nrow(X))
   m <- 2*sum(p*(1-p)) 
  } else m <- 1
  covgbv <- G%x%Phi 
  R <- P-G  
  L <- R%*%solve(G) # Cambie solve() por ginv()
  Vy <-covgbv+R%x%diag(g)
  Z <- diag(t)%x%X
  A <- Z%*%t(Z)
  B <- A%*%solve(A/m+L%x%diag(g))/m # Cambie solve() por ginv()
  y <- apply(allmean,2,scale, scale=F)
  ahat <- B%*%as.vector(as.matrix(y))
  ahat <- matrix(ahat,ncol=t,byrow=FALSE)
  pev <- covgbv%*%ginv(Vy)%*%covgbv
  covgebv <- ginv((diag(t)%x%diag(g))+(G%*%solve(R))%x%Phi)%*%(G%x%Phi)
  #covgebv <- B%*%Vy%*%t(B)
  congebv <- diag(t*g)-ginv(covgbv)%*%pev  
  Gama <- matrix(0,t,t)
  for(i in 1:t){
    for(j in i:t) Gama[i,j] <- mean(diag(covgebv[(1:g)+(i-1)*g,(1:g)+(j-1)*g]))
      #if(i==j) Gama[i,j] <- Gama[i,j]*0.5
  }
  Gama <- t(Gama)+Gama
  diag(Gama) <- 0.5*diag(Gama)
  dimnames(Gama) <- dimnames(covG)
  return(list(cov=covgebv,conf=congebv,pev=pev,ga=Gama,ebv=ahat))
}


