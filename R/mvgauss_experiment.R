#!/usr/bin/env Rscript

# example call
# Rscript R/mvgauss_experiment.R --ms 100 500 2500 --pi1 0.1 --nbhd_size 100 --target_pow 0.6 --rho 0.5 --sidetype one --deptype block --positype uniform --mutype both --nreps 100 --seed 1

devtools::load_all('/Users/destati/Dropbox/001_Personal_Documents/000_Academic_Projects/graphbh/graphbh_dev')

#library(graphbh)
library(here)
source('R/experiment_functions.R')




if (!interactive()){
    suppressPackageStartupMessages(library(argparse))
    
    parser <- ArgumentParser()
    
    parser$add_argument("--ms", type = "integer", nargs = "+", help = "Number of hypotheses")
    parser$add_argument("--pi1", type = "double", help = "Proportion of non-nulls")
    parser$add_argument("--nbhd_size", type = "integer", help = "Size of dependent neighborhood")
    parser$add_argument("--target_pow", type = "double", help = "Target power")
    parser$add_argument("--rho", type = "character", help = "Correlation value, or the string maxneg")
    parser$add_argument("--sidetype", type = "character", help = "Side of the test")
    parser$add_argument("--deptype", type = "character", help = "Whether dependence is block or banded")
    parser$add_argument("--positype", type = "character", help = "Whether signal position is fixed or random")
    parser$add_argument("--mutype", type = "character", help = "Whether signal strength is fixed, random, or both")
    parser$add_argument("--scalability", type = "logical", help = "Whether to do a scalability experiment", default = FALSE)
    
    parser$add_argument("--nreps", type = "integer", help = "Number of repeats")
    parser$add_argument("--seed", type = "integer", help = "Random seed")
    
    parser$add_argument("--signal", type = "double", help = "Signal strength")
    
    parser$add_argument("--notes", type = "character", help = "notes for filename")
    
    
    
    args <- parser$parse_args()
    
    ms <- args$ms
    pi1 <- args$pi1
    nbhd_size <- args$nbhd_size
    target_pow <- args$target_pow
    rho <- ifelse(args$rho == 'maxneg', -Inf, as.numeric(args$rho))
    sidetype <- args$sidetype
    deptype <- args$deptype
    positype <- args$positype
    mutype <- args$mutype
    scalability <- args$scalability
    nreps <- args$nreps
    seed <- args$seed
    signal <- args$signal
    notes <- args$notes

} else {
    ms <- c(10^6)
    pi1 <- 0.1
    nbhd_size <- 100
    target_pow <- 0.6
    rho <- 0.5
    sidetype <- 'two'
    deptype <- 'block'
    positype <- 'uniform'
    mutype <- 'fixed'
    scalability <- TRUE
    nreps <- 10
    seed <- 1
    signal <- 3
    notes <- NULL
}

set.seed(seed)

mstring = paste(ms, collapse="-")

if(mutype == 'both'){
    mutype = c('fixed', 'random')
}


for(mt in mutype){
    filename <- paste0("rawdata/mvgauss",
                       "_ms", mstring, 
                       "_pi1", pi1, 
                       "_nbhd_size", nbhd_size,
                       "_target_pow", target_pow,
                       "_rho", rho,
                       "_sidetype", sidetype,
                       "_deptype", deptype,
                       "_positype", positype,
                       "_mutype", mt,
                       "_scalability", scalability,
                       "_nreps", nreps,
                       "_seed", seed,
                       "_notes", notes,
                       ".RDS")
    calibsig = is.null(signal)
    result = mvgauss_experiment(nreps = nreps,
                                ms = ms,
                                nbhd_size = nbhd_size,
                                alpha = 0.1,
                                calibsig = calibsig,
                                calibeta = TRUE,
                                target_pow = target_pow,
                                pi1 = pi1,
                                rho = rho,
                                sidetype = sidetype,
                                deptype = deptype,
                                positype = positype, 
                                mutype = mt,
                                scalability = scalability,
                                signal = signal)
    saveRDS(result, file = here(filename))
}


