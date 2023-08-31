Readwgt <-
function () 
{
    n <- ntraits
    nom <- wgts[as.logical(wgts[,2]),1]
    ewgt <- wgts[1:n, 3]
    names(ewgt) <- nom
    return(ewgt)
}

