GenRel <- function(markers,sel=c('f2','dh')){  
  X <- markers
  #print(X)
  XXT <- tcrossprod(X)
  m <- ncol(X)
  r <- nrow(X)
  if(sel[1]=='f2'){
    Phi <- XXT*2/m
  } else if(sel[1]=='dh') Phi <- XXT*4/(3*m)-matrix(1,r,r)
  return(Phi)
}



