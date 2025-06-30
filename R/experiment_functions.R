################################################################################
# NOTE: put this at the top of your script, before defining the functions:
################################################################################
library(graphbh)
library(pbapply)
library(parallel)
source('R/utils.R')


################################################################################
# 1) mvgauss_experiment WITH parallelization
################################################################################

mvgauss_experiment <- function(nreps = 1000,
                               ms    = c(20, 100),
                               nbhd_size = 20,
                               alpha  = 0.05,
                               rho    = 0.5,
                               pi1    = 0.2,
                               calibsig   = TRUE,
                               calibeta   = TRUE,
                               target_pow = 0.6,
                               signal     = 3,
                               sidetype   = c('one', 'two'),
                               deptype    = c('block', 'banded'),
                               positype   = c('uniform', 'clustered'),
                               mutype     = c('fixed', 'random'),
                               scalability = FALSE)
{
    # ──────── 1. Match arguments & precompute anything outside the loop ────────
    sidetype  <- match.arg(sidetype)
    deptype   <- match.arg(deptype)
    positype  <- match.arg(positype)
    mutype    <- match.arg(mutype)
    
    Sigmachol <- NA
    if (deptype == 'banded') {
        chols <- lapply(ms, function(m) {
            Sigma_m <- make_banded_cov(m, nbhd_size, rho)
            cholobj  <- Matrix::Cholesky(Sigma_m)
            Matrix::expand(cholobj)$L
        })
        names(chols) <- as.character(ms)
    }
    
    results_all_m <- vector(mode = 'list', length = length(ms))
    names(results_all_m) <- as.character(ms)
    
    thomas_sig <- 6
    lambda0    <- 20
    
    ##############
    # 2. Define the inner helper‐functions (they “close over” sidetype, deptype, etc.)
    ##############
    
    # 2a) generate signal locations
    mvgauss_signal_locs <- function(m, pi1, eta0) {
        if (positype == 'uniform') {
            sample(m, floor(m * pi1))
        } else {
            # clustered via Thomas process
            nonnull_from_thomasPP(m, eta0, lambda0, sigma = thomas_sig)
        }
    }
    
    # 2b) generate p‐values given one draw of Z ∼ MVN
    mvgauss_pvals <- function(m, signal, H1_idx) {
        if (scalability){
            Z <- matrix(rnorm(n = m), nrow = 1)
        }
        else if (deptype == 'block') {
            Z <- rnorm_blk(n = 1, p = m, nbhd_size, rho)
        } 
        else {
            # deptype == 'banded'
            Z <- rmvn(n = 1, mu = rep(0, m), chol.mat = chols[[as.character(m)]])
        }
        
        X <- Z
        if (mutype == 'fixed') {
            mu <- signal
        } else {
            mu <- rexp(length(H1_idx), rate = 1 / signal)
        }
        X[, H1_idx] <- Z[, H1_idx] + mu
        
        x <- as.numeric(X)
        if (sidetype == 'one') {
            pnorm(x, lower.tail = FALSE)
        } else {
            2 * pnorm(abs(x), lower.tail = FALSE)
        }
    }
    
    # 2c) calibrate eta0 so that the Thomas‐process average π₁ ≈ target_pi1
    eta0_calib <- function(m, target_pi1) {
        pi1_diff <- function(eta0) {
            pi1samps <- replicate(1000, {
                H1_idx <- mvgauss_signal_locs(m, pi1, eta0)
                length(H1_idx) / m
            })
            mean(pi1samps) - target_pi1
        }
        find_root(pi1_diff)
    }
    
    # 2d) calibrate “signal” so that BH‐power ≈ target_pow
    signal_calib <- function(m, target_pow, pi1, eta0) {
        power_diff <- function(signal_val) {
            powsamps <- replicate(5000, {
                H1_idx <- mvgauss_signal_locs(m, pi1, eta0)
                pvals   <- mvgauss_pvals(m, signal_val, H1_idx)
                conftab <- results_table(m, method = 'BH', H1_idx, pvals, alpha)
                pow     <- conftab['tp'] / (conftab['tp'] + conftab['fn'])
                if (is.na(pow)) NA else pow
            })
            mean(powsamps, na.rm = TRUE) - target_pow
        }
        
        if (power_diff(0) > 0) {
            return(0)
        } else {
            find_root(power_diff)
        }
    }
    
    
    ##############
    # 3. Main loop over m, but parallelize the “replicate” step inside
    ##############
    
    for (m in ms) {
        m_char <- as.character(m)
        cat("→ Running m =", m_char, "…\n")
        
        # (A) Build adjlist once
        if (deptype == 'block') {
            adjlist <- make_block_adjlist(m, nbhd_size)
        } else {
            adjlist <- make_banded_adjlist(m, nbhd_size)
        }
        
        # (B) calibrate eta0, signal if needed
        eta0 <- m * pi1 / lambda0
        if (calibeta && positype == 'clustered') {
            eta0 <- eta0_calib(m = m, target_pi1 = pi1)
        }
        
        if (calibsig) {
            signal <- signal_calib(m = m,
                                   target_pow = target_pow,
                                   pi1 = pi1,
                                   eta0 = eta0)
        }
        cat("   • Calibrated signal =", round(signal, 2), "\n")
        
        # ──── 3.1) Build a PSOCK cluster and export everything the workers need ────
        ncores <- max(detectCores() - 1, 1)
        #ncores <- 2
        cl     <- makeCluster(ncores)
        
        # Which objects/functions must be copied into each worker’s global env?
        to_export <- c(
            # helper‐functions defined above:
            "mvgauss_signal_locs",
            "mvgauss_pvals",
            "results_table",
            # constants that the helpers close over:
            "m", "pi1", "eta0", "signal", "rho", "scalability",
            # adjacency info, etc.:
            "adjlist", "nbhd_size", "alpha"
        )
        
        clusterExport(cl, varlist = to_export, envir = environment())
        
        # If any of your helpers rely on a loaded package (e.g. graphbh, Matrix),
        # you must load those inside each worker:
        clusterEvalQ(cl, {
            library(graphbh)    # for results_table(…)
            library(Matrix)     # if rmvn() or Cholesky is used (though chols was already computed)
            # (add other library(...) calls if your helpers need them)
            source('R/utils.R')
        })
        
        # ──── 3.2) Now run “nreps replicates” in parallel ─────────────────────────────
        pbapply::pboptions(type = "timer")
        res_list <- pblapply(
            1:nreps,
            function(i) {
                methods     <- c("BH", "eBH", "BYgraph", "BY", "IndBH1", "IndBH2")
                result_list <- vector("list", length(methods))
                names(result_list) <- methods
                
                # 1) generate new H1‐indices + p‐values
                H1_idx <- mvgauss_signal_locs(m = m, pi1 = pi1, eta0 = eta0)
                pvals   <- mvgauss_pvals(m = m, signal = signal, H1_idx = H1_idx)
                
                # 2) run each method
                for (meth in methods) {
                    result_list[[meth]] <- results_table(
                        m, meth, H1_idx, pvals, alpha, adjlist, nbhd_size
                    )
                }
                
                # 3) stack them into a single data.frame (rows = methods)
                do.call(rbind, result_list)
            }#,
            #cl = cl
        )
        
        # ──── 3.3) Shut down the cluster and store results ────────────────────────────
        stopCluster(cl)
        
        # res_list is a list of length nreps, each element a data.frame (#methods × …)
        results_all_m[[m_char]] <- simplify2array(res_list)
    }
    
    return(results_all_m)
}



################################################################################
# 2) adv_experiment WITH parallelization
################################################################################

adv_experiment <- function(nreps     = 1000,
                           ms        = c(20, 100),
                           nbhd_size = 20,
                           alpha     = 0.05,
                           rho       = 0.5,
                           pi1       = 0.2)
{
    results_all_m <- vector(mode = 'list', length = length(ms))
    names(results_all_m) <- as.character(ms)
    
    ##############
    # 2a. Inner helper: uniform signal‐locs
    ##############
    adv_signal_locs <- function(m, pi1) {
        sample(m, floor(m * pi1))
    }
    
    
    ##############
    # 3. Main loop over m; parallelize the replicate inside
    ##############
    for (m in ms) {
        m_char <- as.character(m)
        cat("→ Running adv_experiment, m =", m_char, "…\n")
        
        adjlist <- make_block_adjlist(m, nbhd_size)
        
        # ──── 3.1) Build cluster, export helpers/constants ───────────────────────
        ncores <- max(detectCores() - 1, 1)
        cl     <- makeCluster(ncores)
        
        to_export <- c(
            "adv_signal_locs",
            "results_table",
            "radvblk_p",       # I assume radvblk_p() is defined in utils.R
            # constants closed over by helpers:
            "m", "pi1", "nbhd_size", "alpha",
            "adjlist"
        )
        clusterExport(cl, varlist = to_export, envir = environment())
        
        # load needed packages on each worker
        clusterEvalQ(cl, {
            library(graphbh)    # for results_table()
            # load anything else that radvblk_p(…) or results_table(…) needs
        })
        
        # ──── 3.2) Parallel “replicates” ──────────────────────────────────────────
        pbapply::pboptions(type = "timer")
        res_list <- pblapply(
            1:nreps,
            function(i) {
                methods     <- c("BH", "eBH", "BYgraph", "BY", "IndBH1", "IndBH2")
                result_list <- vector("list", length(methods))
                names(result_list) <- methods
                
                H1_idx <- adv_signal_locs(m = m, pi1 = pi1)
                pvals   <- radvblk_p(n = 1, m = m, nbhd_size = nbhd_size, alpha = alpha)
                
                for (meth in methods) {
                    result_list[[meth]] <- results_table(
                        m, meth, H1_idx, pvals, alpha, adjlist, nbhd_size
                    )
                }
                do.call(rbind, result_list)
            },
            cl = cl
        )
        
        # ──── 3.3) Stop cluster + store ──────────────────────────────────────────
        stopCluster(cl)
        results_all_m[[m_char]] <- res_list
    }
    
    return(results_all_m)
}

################################################################################
# END of parallelized functions
################################################################################
