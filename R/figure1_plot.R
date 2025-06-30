library(ggplot2)
library(here)

file = 'rawdata/wrangled_mvgauss_ms100-500-2500-10000-50000_pi10.1_nbhd_size100_target_pow_rho0.5_sidetypetwo_deptypeblock_positypeuniform_nreps2000_seed1_notesforfig1.RDS'
df4plot = readRDS(file)

# Create annotation data frame for hlines
hline_df <- data.frame(
    metric = c("tprat"),
    yintercept = c(1)  # Different values for each facet
)

# Create relabeler for facet
facet_label <- as_labeller(
    c("fdp" = "FDR", "tprat" = "TP Ratio", "pow" = "Power", 
      "fixed" = "fixed signal strength", "random" = "random signal strength")
)

# Remove eBH (too conservative) and power (using tprat instead)
df4plot = dplyr::filter(df4plot, 
                        method != "BYgraph",
                        method != "BH",
                        method != "IndBH1", 
                        method != "eBH", 
                        metric == "tprat",
                        mutype == 'fixed' )

# Adjust typos
df4plot$method[df4plot$method == "IndBH1"] <- "IndBH"
df4plot$method[df4plot$method == "IndBH2"] <- "IndBH3"


# Plot
plt = ggplot(df4plot, aes(x = m, y = mean, color = method)) + 
    # Add line
    geom_line() + 
    geom_point(size = 0.5) + 
    # Add error bars
    # geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1,
    #               alpha = 0.3, position = position_dodge(width = 0.05)) +
    # Log scale
    scale_x_log10(minor_breaks = c(10^2*(1:10), 10^3*(1:10), 10^4*(1:10))) + 
    # Hlines
    geom_hline(data = hline_df, aes(yintercept = yintercept), linetype = "dashed", color = "black") + 
    # Labels
    labs(x = "m", y = NULL) +
    ylim(0.2, 1) + 
    ylab('TP Ratio') + 
    
    # Theme
    theme_bw() + 
    theme(panel.grid.minor = element_line(color = "grey95", size = 0.3), 
          strip.placement = 'outside') 

# Save the figure, with filename removing wrangled & .RDS
#pltname = gsub("wrangled_|\\.RDS$", "", basename(file))
pltname = 'figure1'
pdffn = paste0('plots/', pltname, '.pdf')
ggsave(pdffn, plot = plt, width = 5, height = 2)