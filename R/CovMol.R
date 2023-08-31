CovMol <-
function(codes,QTLS){
    qtl <- as.matrix(QTLS$QTL)
    link <- as.matrix(QTLS$Link)
    scores <- matrix(0,nrow(codes),ncol(qtl))
    for(i in 1:ncol(qtl)){
    scores[,i] <- codes[,link[!is.na(link[,i]),i]]%*%qtl[!is.na(qtl[,i]),i]
    }
    covar <- cov(scores)
    return(list(S=covar,Scores=scores))
}

