IndCoef <- function(theta,rest=NULL,prest=NULL,Gamma=NULL,G=NULL,evd=FALSE){
  if(!is.null(Gamma)){
    t <- length(theta)
    
    if(!is.null(rest)||!is.null(prest)){
      id.rest <- which(rest==1)
      r <- length(id.rest)
      W <- matrix(0,nrow=t, ncol=r)
      if(!is.null(prest)){ 
        prest <- prest[id.rest]
        C <- rbind(diag(prest[r],r-1),-prest[-r])
        } else C <- diag(r)
      for(j in 1:r){
        W[id.rest[j],j] <- 1
      }
      Q <- (diag(t)-W%*%C%*%solve(t(C)%*%t(W)%*%Gamma%*%W%*%C)%*%t(C)%*%t(W)%*%Gamma)
      } else Q <- diag(t)
    
    if(evd){     
      beta <- theta*abs(svd(Q%*%solve(G)%*%Gamma)$u[,1])
    } else beta <- Q%*%theta  
  } else  beta <- theta
return(as.vector(beta))
}


