library(ggplot2)
library(here)

# List wrangled files
files = list.files(path = here("rawdata"), pattern = "^wrangled_ctrl.*\\.RDS$", full.names = TRUE)

for(file in files){
    
    df4plot = readRDS(file)
    
    # Create annotation data frame for hlines
    hline_df <- data.frame(
        metric = c("fdp", "rejrat"),
        yintercept = c(0.5, 1)  # Different values for each facet
    )
    
    # Create relabeler for facet
    facet_label <- as_labeller(
        c("fdp" = "FDR", "rejrat" = "Rej. Ratio", "pow" = "Power", 
          "gaussian" = "Gaussian", "adversarial" = "Adversarial")
    )
    
    # Remove eBH (too conservative) and power
    df4plot = dplyr::filter(df4plot, 
                            method != "BYgraph",
                            method != "eBH", 
                            metric != "pow",
                            metric != "tprat")
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
        #               alpha = 0.3, position = position_dodge(width = 0.02)) +
        # Log scale
        scale_x_log10(minor_breaks = ) + 
        # Faceting
        facet_grid(metric ~ distype, 
                   scales = "free_y", 
                   labeller = facet_label, 
                   switch="y") + 
        # Hlines
        geom_hline(data = hline_df, aes(yintercept = yintercept), linetype = "dashed", color = "black") +  
        # Labels
        labs(x = "Number of hypotheses (m)", y = NULL) +
        
        # Theme
        theme_bw() + 
        theme(panel.grid.minor = element_line(color = "grey95", size = 0.2), 
              strip.placement = 'outside') 
    
    # Save the figure, with filename removing wrangled & .RDS
    pltname = gsub("wrangled_|\\.RDS$", "", basename(file))
    pdffn = paste0('plots/', pltname, '.pdf')
    ggsave(pdffn, plot = plt, width = 6, height = 3)
}



