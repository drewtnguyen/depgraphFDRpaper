
# Access necessary packages
library(data.table)
library(tidyverse)
library(readxl)
library(snpStats)
library(Matrix)
library(gridExtra)  # for grid.arrange
library(ggplot2)
library(graphbh)
library(igraph)
# devtools::load_all('/Users/destati/Dropbox/001_Personal_Documents/000_Academic_Projects/graphbh/graphbh_dev')

source('R/gwas_utils.R')
source('R/utils.R')


# Load the GWAS dataset 
# from https://pgc.unc.edu/for-researchers/download-results/, 
# leads to https://figshare.com/articles/dataset/scz2014/14672163?file=28570554

scz_data <- fread("gunzip -c realdata/daner_PGC_SCZ52_0513a.hq2.gz", 
                  select = c("CHR", "BP", "SNP", "INFO", "P")) %>%
    .[INFO >= 0.6]
setkey(scz_data, SNP)

# Read in the PLINK data ----- takes several min
g1000_eur_fam_path <- "realdata/g1000_eur/g1000_eur.fam"
g1000_eur_bim_path <- "realdata/g1000_eur/g1000_eur.bim"
g1000_eur_bed_path <- "realdata/g1000_eur/g1000_eur.bed"

ref_genotypes <- read.plink(g1000_eur_bed_path, 
                            g1000_eur_bim_path,
                            g1000_eur_fam_path)$genotypes

# Next, filter out  "eSNP"s from brainvar
brainvar_eqtl_sheets <- read_excel_allsheets("realdata/brainvar_eqtls.xlsx")
scz_data <- scz_data[SNP %in% unique(brainvar_eqtl_sheets$HCP_allEQTLs_FDR0.05$rsID),]

# and keep only those SNPs in both that are included in both scz_data and ref_genotypes
common_snps <- intersect(scz_data$SNP, colnames(ref_genotypes))
scz_data <- scz_data[SNP %in% common_snps]
ref_genotypes <- ref_genotypes[, scz_data$SNP]

# Add columns
scz_data = scz_data %>% 
    mutate(
        CHR = as.factor(CHR),
        neglog10P = -log10(P)
    )

# Now, sample them down to just chromosome 1
scz_data_sub = scz_data %>% 
    filter(CHR == 1)

# Then sort in chromosome order:
scz_data_sub <- scz_data_sub %>%
    arrange(CHR, BP)

# Only include the eSNPs, sorted by location
ref_genotypes_sub <- ref_genotypes[, scz_data_sub$SNP]

# Check the maximum depth 
max_depth = find_depth(scz_data_sub$BP, threshold = 500*10^3)

# free some RAM
rm(brainvar_eqtl_sheets)
#rm(ref_genotypes)

# Now check the LD

LDsnp <- ld(ref_genotypes_sub, stats=c("R.squared"), depth=max_depth)
depmat = (LDsnp > 0.05)
depmat = forceSymmetric(depmat)
depplt = image(depmat, lwd = 0, useRaster = TRUE, col.regions = c("white", "black"), sub = "Dependency Matrix (Chromosome 1)", cex.sub = 0.2)



# Plot
manhatplt = ggplot(scz_data_sub, aes(x = BP, y = neglog10P)) +
    geom_point(alpha = 0.6, color = 'steelblue') +
    labs(
        x = "Genomic Position (Chromosome 1)",
        y = expression(-log[10](p))
    ) +
    theme_minimal() +
    theme(
        legend.position = "none",
        panel.grid.minor.x = element_blank()
    )


# Finally, run the methods
print('Running methods...')
reslist = vector(mode = 'list', length = 2)

for(i in 1:2){
    if(i == 1){
        pvals = scz_data_sub$P
        xlab = "Method (Chromosome 1 p-values)"
        LDsnp <- ld(ref_genotypes_sub, stats=c("R.squared"), depth=max_depth)
    }
    else{
        pvals = scz_data$P
        xlab = "Method (all chromosomes)"
        LDsnp <- ld(ref_genotypes, stats=c("R.squared"), depth=max_depth)
    }
    # first, subset with the BH procedure
    rejBH = BH(alpha = 0.1, pvals)
    pvalsBH = pvals[rejBH]
    # obtain the dependency graph
    depmat = (LDsnp > 0.05)
    depmat = forceSymmetric(depmat)
    
    # blockify
    # ---- 1.  connected–component labels ---------------------------------
    g        <- graph_from_adjacency_matrix(depmat, mode = "undirected")
    memb     <- components(g)$membership          # integer vector length = n
    n        <- length(memb)
    n_comp   <- max(memb)                         # number of components

    ## ---- 2.  one‑hot encode the membership vector -----------------------
    ## This has exactly n non‑zeros, so it is trivially cheap to store.
    Z <- sparseMatrix(i = seq_len(n), j = memb, x = TRUE,
                      dims = c(n, n_comp), repr = "C")   # logical n×n_comp

    ## ---- 3.  same‑component ⇒ 1  (complete graph inside each block) -----
    ## tcrossprod(Z) = Z %*% t(Z) gives an n×n logical matrix whose (i,j)‑entry
    ## is TRUE  ⇔  vertices i and j are in the same component.
    Cliq <- tcrossprod(Z)
    Cliq <- as(Cliq, "lsCMatrix")          # symmetric logical sparse matrix



    depmatBH = Cliq[rejBH, rejBH]
    adjlistBH = adjmat_to_adjlist(depmatBH)
    # depmatBH = depmat[rejBH, rejBH]
    # adjlistBH = adjmat_to_adjlist(depmatBH)
    
    # Now try the other procedures
    RBH = length(rejBH)
    m = length(pvals)
    rejBY = BY(alpha = 0.1, pvals)
    rejIndBH = IndBH(alpha = 0.1*RBH/m, pvalsBH, adjlistBH)
    rejIndBH2 = IndBH_plus(alpha = 0.1*RBH/m, recurse = 1, pvalsBH, adjlistBH)
    
    results = c(
        length(rejBH),
        length(rejBY), 
        length(rejIndBH),
        length(rejIndBH2)
    )
    
    reslist[[i]] = results
    
}


method = c(
    'BH', 
    'BY',
    'IndBH', 
    'IndBH2'
)

resdf1 = data.frame(method = method, results = reslist[[1]], chromosome = 'Chromosome 1 only')
resdfall = data.frame(method = method, results = reslist[[2]], chromosome = 'All Chromosomes')
resdf = rbind(resdf1, resdfall)

rejplt = ggplot(resdf, aes(x = method, y = results, fill = chromosome)) +
    geom_bar(stat = "identity", width = 0.6) +
    scale_fill_manual(
        values = c("Chromosome 1 only" = 'steelblue',
          "All Chromosomes" = 'darkorange2'), 
        name = "SNP Locations"
    ) + 
    theme_minimal() +
    labs(
        x = 'Method',
        y = "Number of Rejections"
    ) +
    theme(
        legend.position = c(0.75, 0.70),          # ⬅️ position inside plot
        legend.background = element_rect(fill = "white", color = "black"),  # optional: make it stand out
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm")  # smaller squares
    )


# Save all to one PDF
pdf("plots/gwas_plot2.pdf", width = 10, height = 3.5)
grid.arrange(manhatplt, depplt, rejplt, ncol = 3)
dev.off()

# # check how long it takes
# start_time <- Sys.time()
# IndBH_plus(alpha = 0.1*RBH/m, recurse = 2, pvalsBH, adjlistBH)
# end_time <- Sys.time()
# end_time - start_time  
# 
# start_time <- Sys.time()
# IndBH(alpha = 0.1*RBH/m, pvalsBH, adjlistBH, block = TRUE)
# end_time <- Sys.time()
# end_time - start_time  








