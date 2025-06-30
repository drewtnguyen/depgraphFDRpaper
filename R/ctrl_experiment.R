#!/usr/bin/env Rscript

# example call

library(graphbh)
library(here)
source('R/experiment_functions.R')


if (!interactive()){
    suppressPackageStartupMessages(library(argparse))
    
    parser <- ArgumentParser()
    
    parser$add_argument("--ms", type = "integer", nargs = "+", help = "Number of hypotheses")
    parser$add_argument("--alpha", type = "double", help = "Level")
    parser$add_argument("--pi1", type = "double", help = "Proportion of non-nulls")
    parser$add_argument("--nbhd_size", type = "integer", help = "Size of dependent neighborhood")

    parser$add_argument("--nreps", type = "integer", help = "Number of repeats")
    parser$add_argument("--seed", type = "integer", help = "Random seed")
    
    
    args <- parser$parse_args()
    
    ms <- args$ms
    alpha <- args$alpha
    pi1 <- args$pi1
    nbhd_size <- args$nbhd_size
    target_pow <- args$target_pow
    rho <- ifelse(args$rho == 'maxneg', -Inf, as.numeric(args$rho))
    sidetype <- args$sidetype
    deptype <- args$deptype
    positype <- args$positype
    mutype <- args$mutype
    nreps <- args$nreps
    seed <- args$seed
    
} else {
    ms <- c(3, 6, 9, 18, 27)
    alpha <- 0.5
    pi1 <- 0
    nbhd_size <- 3
    nreps <- 10000
    seed <- 1
}

set.seed(seed)

mstring = paste(ms, collapse="-")

distype = c('gaussian', 'adversarial')

for(dist in distype){
    filename <- paste0("rawdata/ctrl",
                       "_ms", mstring, 
                       "_nbhd_size", nbhd_size,
                       "_distype", dist,
                       "_nreps", nreps,
                       "_seed", seed,
                       ".RDS")
    
    if(dist == 'gaussian'){
        result = mvgauss_experiment(nreps = nreps,
                                    ms = ms,
                                    nbhd_size = nbhd_size,
                                    alpha = alpha,
                                    calibsig = FALSE,
                                    calibeta = FALSE,
                                    pi1 = 0,
                                    rho = -sqrt(2)/4,
                                    sidetype = 'one',
                                    deptype = 'block',
                                    positype = 'uniform', 
                                    mutype = 'fixed')
    }
    else if(dist == 'adversarial'){
        result = adv_experiment(nreps = nreps,
                                ms = ms,
                                nbhd_size = nbhd_size,
                                alpha = alpha,
                                pi1 = 0)
    }
    saveRDS(result, file = here(filename))
}


