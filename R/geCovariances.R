geCovariances <-
function (data, design = design) 
{
    ncols <- ntraits
    is <- unlist(sapply(1:ncols, function(i) rep(i, ncols - i + 
        1)))
    js <- unlist(sapply(1:ncols, function(i) i:ncols))
    nom <- colnames(data[, -(1:(fcol - 1))])
    adjcovij <- function(i, j, data, design) {
        coupij <- function(i, j, data) {
            data[, i] + data[, j]
        }
        coup12 <- data.frame(data[, 1:(fcol - 1)], coupij(i, 
            j, data))
        cov1 <- Adjcovs2(i, data, design)[2]
        if (i != j) {
            cov2 <- Adjcovs2(j, data, design)[2]
            cov3 <- Adjcovs2(4, coup12, design)[2]
            covar <- (cov3 - cov1 - cov2) * 0.5
        }
        else {
            covar <- cov2 <- cov1
        }
        covij <- c(i - fcol + 1, j - fcol + 1, covar)
        names(covij) <- c("row", "col", "cov")
        return(covij)
    }
    matij <- function(i, j, data) {
        if (i == 1 & j == 1) {
            ncols <- max(data[, 1])
            mat <<- matrix(rep(0, ncols^2), ncol = ncols)
        }
        ind <- which((data[, 1] == i) & (data[, 2] == j))
        mat[i, j] <<- ifelse(i != j, data[ind, 3], data[ind, 
            3] * 0.5)
    }
    covdata <- t(mapply(adjcovij, fcol - 1 + is, fcol - 1 + js, 
        MoreArgs = list(data = data, design = design)))
    mapply(matij, is, js, MoreArgs = list(data = covdata))
    matcov <- mat + t(mat)
    dimnames(matcov) <- list(nom, nom)
    return(matcov)
}

