library(ggplot2)
source('R/utils.R')

set.seed(4)


for(m in c(500, 2500, 10000)){
    pltname = paste0("thomas_rug_", m)
    pdffn = paste0('plots/', pltname, '.pdf')
    
    sizescale = 100/m
    b = 100
    pi0 = 0.9
    lambda0 = 20
    ypos = 0.5
    ypos_label <- 0.55          # vertical position for the label (adjust if needed)
    barendlength = 0.1
    
    
    thomas_sig = 6
    eta0 = (1 - pi0)*m/lambda0
    H1_idx = nonnull_from_thomasPP(m, eta0, lambda0, sigma = thomas_sig)
    
    
    # Put data in a data frame for ggplot
    df <- data.frame(x = H1_idx)
    
    plt = ggplot(df, aes(x = x)) +
        geom_rug(sides = "b", color = "red", size = 1.5*sizescale, alpha = 0.9, length = unit(0.1, "npc")) +
        annotate("segment", x = m/2 - b/2, xend = m/2 + b/2, y = ypos, yend = ypos, size = 0.5) +
        # End caps
        annotate("segment", x = m/2 - b/2, xend = m/2 - b/2, y = ypos - barendlength/2, yend = ypos + barendlength/2, size = 0.5) +
        annotate("segment", x = m/2 + b/2, xend = m/2 + b/2, y = ypos - barendlength/2, yend = ypos + barendlength/2, size = 0.5) +
        scale_x_continuous(limits = c(1, m)) +
        scale_y_continuous(limits = c(0, 1.1)) +
        theme_minimal() + 
        labs(x = "index", y = "") +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), 
              panel.grid.minor = element_line(color = "grey95", size = 0.2), 
              panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_blank()) + 
        geom_label(
            data = data.frame(x = m*0.25, y = ypos_label,
                              lab = paste0("m = ", m)),
            aes(x = x, y = y, label = lab),
            hjust = 1.15,         
            vjust = 0.5,
            size  = 4,
            label.padding = unit(0.15, "lines"),
            label.r = unit(0.15, "lines")  
        ) 
    print(length(H1_idx)/m)
    
    
    ggsave(pdffn, plot = plt, width = 6, height = 1) 
}


