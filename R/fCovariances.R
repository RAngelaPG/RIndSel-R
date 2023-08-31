fCovariances <-function (data) 
{
    Cov <- var(data,na.rm=TRUE)
    nom <- colnames(data)
    dimnames(Cov) <- list(nom, nom)
    return(Cov)
}

