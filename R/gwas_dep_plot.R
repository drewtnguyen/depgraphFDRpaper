# Access necessary packages
library(data.table)
library(tidyverse)
library(readxl)
library(snpStats)
library(Matrix)
library(gridExtra)  # for grid.arrange
library(ggplot2)
library(igraph)
# Load the library
library(depgraphFDR)

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


# Make manhattan plot
# based on https://r-graph-gallery.com/101_Manhattan_plot.html
don <- scz_data %>% 
    # Compute chromosome size
    group_by(CHR) %>% 
    summarise(chr_len=max(BP)) %>% 
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
    # Add this info to the initial dataset
    left_join(scz_data, ., by=c("CHR"="CHR")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHR, BP) %>%
    mutate( BPcum=BP+tot) %>% 
    mutate(CHRfac = factor(CHR, levels = sort(unique(CHR))))

axisdf = don %>%
    group_by(CHR) %>%
    summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

manhatplt = ggplot(don, aes(x=BPcum, y=-log10(P))) +
    
    # Show all points
    geom_point( aes(color=CHRfac), alpha=0.5, size=0.3) +
    scale_color_manual(values = rep(c("coral", "skyblue"), 11 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
    
    # Custom the theme:
    theme_bw() +
    theme( 
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
    ) + 
    labs(x = "Chromosomal position", y = expression(-log[10](p)))

# Make dependency matrix plots

chrs = c(1, 2)
deppltlist = list()
for(chr in chrs){
    scz_data_sub = scz_data %>% 
        filter(CHR == chr)
    
    # Then sort in chromosome order:
    scz_data_sub <- scz_data_sub %>%
        arrange(CHR, BP)
    
    # Only include the eSNPs, sorted by location
    ref_genotypes_sub <- ref_genotypes[, scz_data_sub$SNP]
    
    # Check the maximum depth 
    max_depth = find_depth(scz_data_sub$BP, threshold = 500*10^3)
    
    # Now make check the LD figures
    
    LDsnp <- ld(ref_genotypes_sub, stats=c("R.squared"), depth=max_depth)
    depmat <- (LDsnp > 0.2)
    depmat <- forceSymmetric(depmat)
    diag(depmat) <- TRUE
    chrstring = paste('Chromosome', chr)
    depplt <- Matrix::image(depmat, lwd = 0, useRaster = TRUE, col.regions = c("white", "black"), sub = paste0("Dep. Matrix (", chrstring, ")"), cex.sub = 0.2)
    deppltlist[[chrstring]] = depplt
}

# Finally, run the methods
alpha = 0.01
print('Running methods...')

pvals = scz_data$P
xlab = "Method (all chromosomes)"
LDsnp <- ld(ref_genotypes, stats=c("R.squared"), depth=max_depth)
# first, subset with the BH procedure
rejBH = BH(alpha = alpha, pvals)
pvalsBH = pvals[rejBH]
# obtain the dependency graph
depmat = (LDsnp > 0.2)
depmat = forceSymmetric(depmat)

# blockify---i.e., enforce block dependence so that FDR control is guaranteed

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


# Now try the other procedures
RBH = length(rejBH)
m = length(pvals)
rejBY = BY(alpha = alpha, pvals)

rejIndBH = IndBH(alpha = alpha*RBH/m, pvalsBH, adjlistBH, block = TRUE)
rejIndBH3 = IndBH_plus(alpha = alpha*RBH/m, recurse = 2, pvalsBH, adjlistBH, block = TRUE)

results = c(
    length(rejBH),
    length(rejBY), 
    length(rejIndBH),
    length(rejIndBH3)
)

method = c(
    'BH', 
    'BY',
    'IndBH1', 
    'IndBH3'
)

rejplt = ggplot(mapping = aes(x = method, y = results)) +
    geom_bar(stat = "identity", width = 0.6, fill = 'steelblue') +
    theme_minimal() +
    labs(
        x = bquote("Method (" * alpha == .(alpha) * ")"),
        y = "Number of Rejections"
    ) +
    theme(
        legend.position = c(0.75, 0.75),          # ⬅️ position inside plot
        legend.background = element_rect(fill = "white", color = "black"),  # optional: make it stand out
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm")  # smaller squares
    )


# Save all to one PDF
pdf("plots/gwas_plot.pdf", width = 8.5, height = 5)
grid.arrange(manhatplt, deppltlist[['Chromosome 1']], deppltlist[['Chromosome 2']], rejplt, layout_matrix = rbind(c(1,1,1), c(2,3,4)), heights = c(1,1.5))
dev.off()

