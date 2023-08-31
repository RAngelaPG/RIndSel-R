Means <-function (data, design=design, softR=softR, one.env=one.env, block.ex=block.ex){
    am <- sapply(fcol:ncol(data), adjmeans, data, design, softR, one.env, block.ex)
    am <- as.data.frame(am)
    colnames(am) <- colnames(data[, -(1:(fcol - 1))])
    return(am)
}

