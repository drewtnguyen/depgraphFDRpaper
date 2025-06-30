# Discrete gaussian
dnorm_discr <- function(i, mu, sigma = 1){
    # Omit mu from the denominator due to location invariance
    dnorm(i, mean = mu, sd = sigma) / sum(dnorm(-500:500, sd = sigma))
}


# Generate locations of non-null hypotheses from the thomas point process, # by first generating locations of cluster centers
# and then of the daughter points. 
nonnull_from_thomasPP <- function(m, eta0, lambda0, sigma){
    # Generate locations of cluster centers
    N_clu = rpois(1, lambda = eta0)
    if(N_clu == 0){
        return(integer(0))
    }
    loc_clu = sample(1:m, N_clu, replace = TRUE)
    # Generate locations of daughter points
    loc_dau_lst = vector(mode = 'list', length = length(loc_clu))
    for(k in 1:length(loc_clu)){
        loc = loc_clu[[k]]
        N_dau = rpois(1, lambda = lambda0)
        # Make a large set of integers to sample from
        ints = floor(loc-20*sigma):ceiling(loc+20*sigma)
        # Define probabilities
        prob = dnorm_discr(ints, mu = loc, sigma = sigma)
        loc_dau_lst[[k]] = sample(ints, size = N_dau, prob = prob, replace = TRUE)
    }
    # Join all the results
    loc_dau = do.call(c, loc_dau_lst)
    # Return those unique elements that are in 1:m
    intersect(1:m, unique(loc_dau))
}

# Returns a named vector
results_table <- function(m, method = "BH", H1, 
                          pvals, alpha, adjlist, nbhd_size){
    start <- Sys.time()
    if(method == "BH"){
        rej = BH(alpha, pvals)
    }
    else if(method == "BY"){
        rej = BY(alpha, pvals)
    }
    else if(method == "BYgraph"){
        alphac = corrected_alpha(alpha, nbhd_size, m)
        rej = BH(alphac, pvals)
    }
    else if(method == "eBH"){
        rej = eBH(alpha, pvals)
    }
    else if(method == 'IndBH1'){
        if(length(BH(alpha, pvals)) == 0){
            rej = integer(0)
        }
        else{
            rej = suppressMessages({IndBH(alpha = alpha,
                        pvals = pvals, adjlist = adjlist)})
        }
    }
    else if(method == 'IndBH2'){
        if(length(BH(alpha, pvals)) == 0){
            rej = integer(0)
        }
        else{
            rej = suppressMessages({IndBH_plus(recurse = 2, alpha = alpha,
                        pvals = pvals, adjlist = adjlist)})
        }
    }
    runtime <- Sys.time() - start
    browser()
    H0 = setdiff(1:m, H1)
    nonrej = setdiff(1:m, rej)
    tp = length(intersect(rej, H1))
    fp = length(intersect(rej, H0))
    fn = length(intersect(nonrej, H1))
    tn = length(intersect(nonrej, H0))
    out = c(tp, fp, fn, tn, runtime)
    names(out) = c('tp', 'fp', 'fn', 'tn', 'runtime')
    return(out)
}

# Modification of the code from Karl Broman (https://github.com/kbroman/broman)
rmvn <- function(n, mu = 0, Sigma = NA, chol.mat = NA) {
    if(any(is.na(Sigma)) & any(is.na(chol.mat))){
        stop("Must Input one of Sigma or chol.mat")
    }

    p <- length(mu)
    if (any(is.na(match(dim(Sigma), p)))) {
        stop("Dimension problem!")
    }
    if(any(is.na(chol.mat))){
        D <- chol(Sigma)
    }else{
        D <- chol.mat
    }
    matrix(rnorm(n * p), ncol = p) %*% D + rep(mu, rep(n, p))
}

make_block_adjlist <- function(m, nbhd_size){
    L = vector(mode = 'list', length = m)
    for (i in 1:(m / nbhd_size)) {
        block <- (nbhd_size * (i - 1) + 1):(nbhd_size * i)
        for (j in 1:nbhd_size) {
            L[[(i-1) * nbhd_size + j]] <- block
        }
    }
    return(L)
}

make_banded_adjlist <- function(m, nbhd_size){
    L = vector(mode = 'list', length = m)
    for (i in 1:m) {
        lower = max(i - floor(nbhd_size/2), 1)
        upper = min(i + floor(nbhd_size/2), m)
        L[[i]] <- lower:upper
    }
    return(L)
}

make_banded_cov <- function(m, nbhd_size, rho) {
    half_band  <- floor(nbhd_size / 2)        
    diag_offsets <- 0:half_band 
    # for each offset from the diagonal, build the diagonal
    diags  <- lapply(diag_offsets, function(d) {
        len <- m - d        # length shrinks near the ends
        rep(rho^d, len) # if d = 0, this is the diagonal.
    })
    Matrix::bandSparse(n = m, k = diag_offsets, diagonals = diags, symmetric = TRUE)
}


adjmat_to_adjlist <- function(M){
    # M must be a square, symmetric matrix
    p = nrow(M)
    L <- vector(mode = 'list', length = p)
    for(i in 1:p){
        L[[i]] = which(as.logical(M[,i]))
    }
    return(L)
}



# returns a matrix
rnorm_equi <- function(n, p, rho, chol.mat = NA){
    if(rho < 0){
        Sigma = (1 - rho) * diag(p) + rho * 1
        X = rmvn(n, mu = numeric(p), Sigma = Sigma, chol.mat = chol.mat)
    }
    else{
        Z = rnorm(n, mean = 0, sd = sqrt(rho))
        E = matrix(rnorm(p * n, mean = 0, sd = sqrt(1 - rho)), n)
        X = Z + E
        X
    }
}


# returns a matrix
rnorm_blk <- function(n, p, nbhd_size, rho, chol.mat = NA){
    if(p %% nbhd_size != 0) stop('exact multiples only, p = nbhd_size * k')
    k = p/nbhd_size
    l = vector(mode = 'list', length = k)
    for(blk in 1:k){
        if(rho == -Inf){
            Xblk = rnorm_equi_maxneg(n, nbhd_size)
        }
        else{
            Xblk = rnorm_equi(n, nbhd_size, rho, chol.mat)
        }
        l[[blk]] = Xblk
    }
    # cbind correlated blocks together
    X = do.call(cbind, l)
    X
}


corrected_alpha <- function(alpha, B, m){
    Lb = sum(1/(1:B))
    alpha_c = (1/Lb) * (m/B) * (1 - (1 - alpha)^(B/m))
    # approximately (1/Lb) * alpha.
    alpha_c
}

BH <- function(alpha, pvals, m = NULL){
    if(is.null(m)) {
        m = length(pvals)
    }
    R = max(which(sort(pvals) <= alpha*(1:length(pvals))/m), -Inf)
    # R returns -Inf
    which(pvals <= alpha * R/m)
}

BY <- function(alpha, pvals){
    m = length(pvals)
    newalpha = alpha/sum(1/(1:m))
    BH(newalpha, pvals)
}

eBH <- function(alpha, pvals){
    newpvals = 2*sqrt(pvals)
    newalpha = sqrt(2*alpha)
    BH(newalpha, newpvals)
}

# Root finding utility
find_root <- function(f, lower = 0, upper = 10, max_upper = 50000) {
    while(upper < max_upper){
        root <- tryCatch(
            uniroot(f, c(lower, upper))$root,
            error = function(e) NA
        )
        
        if (!is.na(root)) {
            return(root)
        }
        
        upper <- upper * 2
    }
    stop('no root found')
    return(NA)
}

# Generate adversarial distribution for one block
radvblk_p_oneblock <- function(n = 1, m, nbhd_size, alpha){
    # return a vector of length nbhd_size
    l = vector(mode = 'list', length = n)
    # generate vector of probabilities for k
    pk = alpha*(1/(1:nbhd_size) * nbhd_size/m)
    if(sum(pk) >= 1) stop("alpha L_m not < 1")
    pk = c(1 - sum(pk), pk)
    # use a multinomial to generate k
    k = {
        x = rmultinom(n, 1, pk)
        apply(x, 2, function(v) which(v == 1))
    }
    k = k - 1 # subtract 1 due to 1-indexing
    # Generate p-values
    # put exactly k of them
    # in intvl I_k = [(k - 1)*alpha/m, k alpha/m]
    # and the rest in [(nbhd_size*alpha/m), 1]
    # if K = 0, they all end up there.
    for(i in 1:n){
        p1 = runif(k[i], min = (k[i] - 1)*alpha/m, max = k[i]*alpha/m)
        p2 = runif(nbhd_size - k[i], min = nbhd_size*alpha/m, max = 1)
        #p2 = runif(nbhd_size - k[i], min = alpha, max = 1)
        p = c(p1, p2)
        # Shuffle the pvalues
        #browser()
        p = sample(p)
        l[[i]] = p
    }
    do.call(rbind, l)
}

radvblk_p <- function(n, m, nbhd_size, alpha){
    if(m %% nbhd_size != 0) stop('exact multiples only, m = nbhd_size * k')
    k = m/nbhd_size
    l = vector(mode = 'list', length = k)
    for(blk in 1:k){
        Xblk = radvblk_p_oneblock(n = n, m, nbhd_size, alpha)
        l[[blk]] = as.data.frame(Xblk)
    }
    # cbind correlated blocks together
    X = do.call(cbind, l)
    if(n == 1){
        X = as.numeric(X)
    }
    X
}

