library(ggplot2)
library(here)

# List wrangled files
files = list.files(path = here("rawdata"), pattern = "^wrangled_mvgauss.*\\.RDS$", full.names = TRUE)

for(file in files){
    
    df4plot = readRDS(file)
    
    # Create annotation data frame for hlines
    hline_df <- data.frame(
        metric = c("fdp", "tprat"),
        yintercept = c(0.1, 1)  # Different values for each facet
    )
    
    # Create relabeler for facet
    facet_label <- as_labeller(
        c("fdp" = "FDR", "tprat" = "TP Ratio", "pow" = "Power", 
          "fixed" = "fixed signal strength", "random" = "random signal strength",
          "runtimes" = "Runtime")
    )
    
    # Remove eBH (too conservative) and power (using tprat instead)
    df4plot = dplyr::filter(df4plot, 
                            # method != "BYgraph",
                            # method != "eBH", 
                            metric != "pow",
                            metric != "runtimes")
    
    # Adjust typos
    df4plot$method[df4plot$method == "IndBH1"] <- "IndBH1"
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
        # Faceting
        facet_grid(metric ~ mutype, 
                   scales = "free_y", 
                   labeller = facet_label, 
                   switch="y") + 
        # Hlines
        geom_hline(data = hline_df, aes(yintercept = yintercept), linetype = "dashed", color = "black") +  
        # Labels
        labs(x = "m", y = NULL) +
        
        # Theme
        theme_bw() + 
        theme(panel.grid.minor = element_line(color = "grey95", size = 0.2), 
              strip.placement = 'outside') 
    
  # Save the figure, with filename removing wrangled & .RDS
    pltname = gsub("wrangled_|\\.RDS$", "", basename(file))
    pdffn = paste0('plots/', pltname, '.pdf')
    ggsave(pdffn, plot = plt, width = 6, height = 3)
}