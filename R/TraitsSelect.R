TraitsSelect <-
function (data, ind) 
{
    traits <- data[, -(1:(fcol - 1))][, as.logical(ind)]
    cat("Traits to be used", "\n")
    cat(colnames(traits), "\n")
    traits <- cbind(data[, (1:(fcol - 1))], traits)
    return(traits)
}

