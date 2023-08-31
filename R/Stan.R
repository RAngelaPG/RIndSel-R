Stan <-
function (data, design = "lattice") 
{
    fstd <- function(i) {
        stdres(lm(data[, i] ~ 1))
    }
    std <- sapply(fcol:ncol(data), fstd)
    colnames(std) <- paste("V", 1:ncol(std), sep = "")
    std <- data.frame(data[, 1:(fcol - 1)], std)
    return(std)
}

