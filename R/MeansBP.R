MeansBP <-function (data, design=design, softR=softR, one.env=one.env, block.ex=block.ex){
    amb <- sapply(fcol:ncol(data), AdjBlup, data, design, softR, one.env, block.ex)
    amb <- as.data.frame(amb)
    colnames(amb) <- colnames(data[, -(1:(fcol - 1))])
    return(amb)
}
