library(tidyverse)
library(here)
source('R/wrangle_functions.R')


# List rawdata files
files = list.files(path = here("rawdata"), pattern = "^ctrl_.*\\.RDS$", full.names = TRUE)

# Extract filenames without path
filenames <- basename(files)

file_keys <- gsub("_distype(gaussian|adversarial)", "", filenames)

# Remove ".RDS" extension
clean_basenames <- gsub("\\.RDS$", "", file_keys)

# Create a named list to store pairs
file_pairs <- list()

# Loop through unique keys and pair files
for (basename in unique(clean_basenames)) {
    matched_files <- files[clean_basenames == basename]
    
    # Ensure there are exactly two files in the pair
    if (length(matched_files) == 2) {
        file_pairs[[basename]] <- matched_files  # Assign using clean basename
    }
}

basenames = names(file_pairs)

# Process dataframe for all file pairs

for (basename in basenames){
    pair = file_pairs[[basename]]
    distypes = c('adversarial', 'gaussian')
    names(pair) = distypes
    # Initialize dataframes
    dflist = list()

    for(distype in distypes){
        file = pair[distype]
        # Load `result` object
        result = readRDS(file)
        
        # Main wrangling logic
        ms = names(result)
        nreps = dim(result[[ms[1]]])[3]
        
        
        # Process into a list for each m 
        tp_ratio_list <- setNames(vector("list", length(ms)), ms)
        rej_ratio_list <- setNames(vector("list", length(ms)), ms)
        fdp_list <- setNames(vector("list", length(ms)), ms)
        pow_list <- setNames(vector("list", length(ms)), ms)
        
        
        for (m in ms) {
            tp_ratios <- sapply(1:nreps, function(n) {
                confarr <- result[[m]][,,n]
                tp_ratio_wrt_BH(confarr)
            })
            
            tp_ratio_list[[m]] <- list(
                mean = rowMeans(tp_ratios, na.rm = TRUE),
                sd = apply(tp_ratios, 1, function(x) sd(x, na.rm = TRUE))
            )
            
            rej_ratios <- sapply(1:nreps, function(n) {
                confarr <- result[[m]][,,n]
                rej_ratio_wrt_BH(confarr)
            })
            
            rej_ratio_list[[m]] <- list(
                mean = rowMeans(rej_ratios, na.rm = TRUE),
                sd = apply(rej_ratios, 1, function(x) sd(x, na.rm = TRUE))
            )
            
            fdps <- sapply(1:nreps, function(n) {
                confarr <- result[[m]][,,n]
                fdp(confarr)
            })
            
            fdp_list[[m]] <- list(
                mean = rowMeans(fdps),
                sd = apply(fdps, 1, sd)
            )
            
            pows <- sapply(1:nreps, function(n) {
                confarr <- result[[m]][,,n]
                pow(confarr)
            })
            
            pow_list[[m]] <- list(
                mean = rowMeans(pows, na.rm = TRUE),
                sd = apply(pows, 1, function(x) sd(x, na.rm = TRUE))
            )
            
        }
        
        # Process into a dataframe with columns 
        # `method`, `m`, `tp_ratio`, `std_dev`
        df <- do.call(rbind, lapply(names(tp_ratio_list), function(m) {
            methods <- names(tp_ratio_list[[m]]$mean)
            data.frame(
                m = as.numeric(m),
                method = methods,
                tprat_mean = tp_ratio_list[[m]]$mean[methods],
                tprat_sd = tp_ratio_list[[m]]$sd[methods],
                rejrat_mean = rej_ratio_list[[m]]$mean[methods],
                rejrat_sd = rej_ratio_list[[m]]$sd[methods],
                fdp_mean = fdp_list[[m]]$mean[methods],
                fdp_sd = fdp_list[[m]]$sd[methods],
                pow_mean = pow_list[[m]]$mean[methods],
                pow_sd = pow_list[[m]]$sd[methods],
                row.names = NULL
            )
        }))
        
        # transform chosen columns into `metric`, `mean`, `sd` columns
        df_long <- df %>%
            tidyr::pivot_longer(cols = c(tprat_mean, rejrat_mean, fdp_mean, pow_mean,
                                         tprat_sd, rejrat_sd, fdp_sd, pow_sd), 
                                names_to = c("metric", ".value"), 
                                names_sep = "_") 
        df_long$distype = distype
        dflist[[distype]] = df_long
    }
    
    df4plot = do.call(rbind, dflist)
    saveRDS(df4plot, file = paste0("rawdata/wrangled_", basename, ".RDS"))
}



