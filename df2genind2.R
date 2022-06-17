df2genind2 <- function (X, sep = NULL, ncode = NULL, ind.names = NULL, loc.names = NULL, pop = NULL, missing = NA, ploidy = 2, type = c("codom", "PA")) {
    if (is.data.frame(X)) 
        X <- as.matrix(X)
    if (!inherits(X, "matrix")) 
        stop("X is not a matrix")
    res <- list()
    type <- match.arg(type)
    n <- nrow(X)
    nloc <- ncol(X)
    ploidy <- as.integer(ploidy)
    if (ploidy < 1L) 
        stop("ploidy cannot be less than 1")
    if (is.null(ind.names)) {
        ind.names <- rownames(X)
    }
    if (is.null(loc.names)) {
        loc.names <- colnames(X)
    }
    if (!is.null(pop)) {
        if (length(pop) != n) 
            stop("length of factor pop differs from nrow(X)")
        pop <- as.factor(pop)
    }
    if (toupper(type) == "PA") {
        mode(X) <- "numeric"
        rownames(X) <- ind.names
        colnames(X) <- loc.names
        temp <- apply(X, 2, function(c) all(is.na(c)))
        if (any(temp)) {
            X <- X[, !temp]
            warning("entirely non-type marker(s) deleted")
        }
        temp <- apply(X, 1, function(r) all(is.na(r)))
        if (any(temp)) {
            X <- X[!temp, , drop = FALSE]
            pop <- pop[!temp]
            warning("entirely non-type individual(s) deleted")
        }
        temp <- apply(X, 2, function(loc) length(unique(loc[!is.na(loc)])) == 
            1)
        if (any(temp)) {
            X <- X[, !temp, drop = FALSE]
            warning("non-polymorphic marker(s) deleted")
        }
        prevcall <- match.call()
        res <- genind(tab = X, pop = pop, prevcall = prevcall, 
            ploidy = ploidy, type = "PA")
        return(res)
    }
    mode(X) <- "character"
    if (is.null(sep)) {
        if (!is.null(ncode)) {
            temp <- nchar(X[!is.na(X)])
            if (ncode < max(temp)) 
                stop("some character strings exceed the provided ncode.")
        }
        if (is.null(ncode)) {
            temp <- nchar(X[!is.na(X)])
            ncode <- max(temp)
        }
        if ((ncode%%ploidy) > 0) 
            stop(paste(ploidy, "alleles cannot be coded by a total of", 
                ncode, "characters", sep = " "))
    }
    tempX <- X
    if (!is.null(sep)) 
        tempX <- gsub(sep, "", X) # 1.7 sec
    tempX <- gsub("^0*$", NA, tempX) # 1.7 sec
    tempX <- gsub("(NA)+", NA, tempX) # 1.7 sec
    temp <- apply(tempX, 2, function(c) all(is.na(c)))
    if (any(temp)) {
        X <- X[, !temp, drop = FALSE]
        tempX <- tempX[, !temp, drop = FALSE]
        loc.names <- loc.names[!temp]
        nloc <- ncol(X)
        warning("entirely non-type marker(s) deleted")
    }
    temp <- apply(tempX, 1, function(r) all(is.na(r)))
    if (any(temp)) {
        X <- X[!temp, , drop = FALSE]
        tempX <- tempX[!temp, , drop = FALSE]
        pop <- pop[!temp]
        ind.names <- ind.names[!temp]
        n <- nrow(X)
        warning("entirely non-type individual(s) deleted")
    }
    n <- nrow(X)
    X[is.na(tempX)] <- NA
    fillWithZero <- function(M, targetN) {
        naIdx <- is.na(M)
        keepCheck <- any(nchar(M) < targetN)
        while (keepCheck) {
            mat0 <- matrix("", ncol = ncol(M), nrow = nrow(M))
            mat0[nchar(M) < targetN] <- "0"
            M <- matrix(paste(mat0, M, sep = ""), nrow = nrow(mat0))
            keepCheck <- any(nchar(M) < targetN)
        }
        M[naIdx] <- NA
        return(M)
    }
    if (is.null(sep) | ploidy == as.integer(1)) {
        X <- fillWithZero(X, targetN = ncode)
        splitX <- list()
        for (i in 1:ploidy) {
            splitX[[i]] <- substr(X, 1, ncode/ploidy)
            X <- sub(paste("^.{", ncode/ploidy, "}", sep = ""), 
                "", X)
        }
    }
    if (!is.null(sep)) {
        if (ploidy > 1) {
            
            
            # system.time(temp <- t(as.matrix(as.data.frame( strsplit(X, sep)))))
            
            # changed to
            na_fill <- function(x) {
              if(is.na(x[1])) return(c(NA,NA))
              return(x)
            }
                          
            temp <- ldply(strsplit(X, sep),na_fill)
            
            splitX <- list()
            for (i in 1:ncol(temp)) {
                splitX[[i]] <- matrix(temp[, i], nrow = n)
            }
        }
        else {
            splitX <- list()
            splitX[[1]] <- X
        }
        temp <- unlist(splitX)
        temp <- temp[!is.na(temp)]
        ncode <- max(nchar(temp)) * ploidy
        splitX <- lapply(splitX, function(Y) fillWithZero(Y,targetN = ncode/ploidy))
    }
    loc.all <- list()
    for (i in 1:nloc) {
        temp <- unlist(lapply(splitX, function(e) e[, i]))
        loc.all[[i]] <- sort(unique(temp[!is.na(temp)]))
    }
    names(loc.all) <- loc.names
    temp <- lapply(1:nloc, function(i) matrix(0, nrow = n, ncol = length(loc.all[[i]]), 
        dimnames = list(NULL, loc.all[[i]])))
    names(temp) <- loc.names
    findall <- function(cha, loc.all) {
        if (is.na(cha)) return(NULL)
        return(which(cha == loc.all))
    }
    for (k in 1:ploidy) {
        for (i in 1:n) {
            for (j in 1:nloc) {
                allIdx <- findall(splitX[[k]][i, j], loc.all[[j]])
                temp[[j]][i, allIdx] <- temp[[j]][i, allIdx] + 1
                if (is.null(allIdx)) {
                  temp[[j]][i, ] <- NA
                }
            }
        }
    }
    nall <- unlist(lapply(temp, ncol))
    loc.rep <- rep(names(nall), nall)
    col.lab <- paste(loc.rep, unlist(loc.all, use.names = FALSE), sep = ".")
    mat <- matrix(unlist(temp), nrow = nrow(temp[[1]]))
    mat <- mat/ploidy
    colnames(mat) <- col.lab
    rownames(mat) <- ind.names
    if (!is.na(missing)) {
        if (missing == 0) {
            mat[is.na(mat)] <- 0
        }
        if (toupper(missing) == "MEAN") {
            moy <- apply(mat, 2, function(c) mean(c, na.rm = TRUE))
            for (j in 1:ncol(mat)) {
                mat[, j][is.na(mat[, j])] <- moy[j]
            }
        }
    }
    prevcall <- match.call()
    res <- genind(tab = mat, pop = pop, prevcall = prevcall, 
        ploidy = ploidy, type = type)
    return(res)
}


