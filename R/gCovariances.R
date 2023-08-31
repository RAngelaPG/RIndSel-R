gCovariances <-function (datos, design, softR, one.env, block.ex){
    ncols <- ntraits
    is <- unlist(sapply(1:ncols, function(i) rep(i, ncols - i + 1)))
    js <- unlist(sapply(1:ncols, function(i) i:ncols))
    nom <- colnames(datos[, -(1:(fcol - 1))])
##################################################################################################################
	adjcovij <- function(i, j, datos, design, softR, one.env, block.ex) {
        coupij <- function(i, j, datos) {
            datos[, i] + datos[, j]
        }
        coup12 <- data.frame(datos[, 1:(fcol - 1)], coupij(i,j, datos))
        cov1 <- Adjcovs(i, data=datos, design=design, softR=softR, one.env=one.env, block.ex=block.ex)["Entry"]
        if (i != j) {
            cov2 <- Adjcovs(j, data=datos, design=design, softR=softR, one.env=one.env, block.ex=block.ex)["Entry"]
            cov3 <- Adjcovs(fcol, data=coup12, design=design, softR=softR, one.env=one.env, block.ex=block.ex)["Entry"]
            covar <- (cov3 - cov1 - cov2) * 0.5
        }
        else {
            covar <- cov2 <- cov1
        }
        covij <- c(i - fcol + 1, j - fcol + 1, covar)
        names(covij) <- c("row", "col", "cov")
        return(covij)
    }
##################################################################################################################	
    matij <- function(i, j, datos) {
        if (i == 1 & j == 1) {
            ncols <- max(datos[, 1])
            mat <<- matrix(rep(0, ncols^2), ncol = ncols)
            R <<- matrix(rep(0, ncols^2), ncol = ncols)
        }
        ind <- which((datos[, 1] == i) & (datos[, 2] == j))
        mat[i, j] <<- ifelse(i != j, datos[ind, 3], datos[ind, 3] * 0.5)
    }
##################################################################################################################
    covdata <- t(mapply(adjcovij, fcol - 1 + is, fcol - 1 + js, MoreArgs = list(datos=datos, design=design, softR=softR, one.env=one.env, block.ex=block.ex)))
	mapply(matij, is, js, MoreArgs = list(datos = covdata))
    matcov <- mat + t(mat)
    dimnames(matcov) <- list(nom, nom)
    rm(mat, pos = ".GlobalEnv")
    
	return(matcov)
}

