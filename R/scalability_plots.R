library(ggplot2)
library(here)

file = 'rawdata/wrangled_mvgauss_ms100000_pi10.1_nbhd_size100_target_pow0_rho0_sidetypeone_deptypeblock_positypeuniform_scalabilityTRUE_nreps5_seed1_notes.RDS'
df4plot = readRDS(file)


# Create relabeler for facet
facet_label <- as_labeller(
    c("fdp" = "FDR", "tprat" = "TP Ratio", "pow" = "Power", 
      "fixed" = "fixed signal strength", "random" = "random signal strength",
      "runtimes" = "Runtime")
)

# Remove eBH (too conservative) and power (using tprat instead)
df4plot = dplyr::filter(df4plot, 
                        method != "BYgraph",
                        method != "eBH", 
                        metric == "runtimes",
                        mutype == 'fixed' )

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
    # Labels
    labs(x = "m", y = NULL) +
    ylim(0.2, 1) + 
    ylab('Runtime') + 
    
    # Theme
    theme_bw() + 
    theme(panel.grid.minor = element_line(color = "grey95", size = 0.3), 
          strip.placement = 'outside') 

# Save the figure, with filename removing wrangled & .RDS
#pltname = gsub("wrangled_|\\.RDS$", "", basename(file))
pltname = 'scalability'
pdffn = paste0('plots/', pltname, '.pdf')
ggsave(pdffn, plot = plt, width = 5, height = 2)