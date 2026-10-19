##----------------------------------------------------------------------------------
## Function to perform multiple fishers enrichment test on a wide format data frame
##----------------------------------------------------------------------------------

## libraries
  library(tidyverse)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(GenomicRanges)
  
  set.seed(1234)
  
  source("/Users/ana207/Desktop/R_Projects/Scenic-results/functions/bed_overlaps.R")
  source("/Users/ana207/Desktop/R_Projects/Scenic-results/functions/pairwise_fishers.R")
  

## Arguments:

#      1. wide_df: A wide format df--usually output from enrichment_wide_df function
#                  Should have a column named "marker_peak" a binary column that says if it is a marker or not

## function:

  enrichment_fishers <- function(wide_df= df_1){
    
## looping through each variable name 
  var_names <- setdiff(colnames(wide_df), c("peak", "marker_peak", "NA"))

# performing 
  loop_results <- list()
  
  for (x in var_names) {
    loop_results[[x]] <- pairwise_fishers(df_data = wide_df,
                                              col_1 = "marker_peak",
                                              col_2 = x,
                                              col_name = x)
    
  }

# binding rows and FDR
  fishers_res <- bind_rows(loop_results) %>%
    mutate(adj_p_val = p.adjust(p_val, method = "BH"))

  return(fishers_res)
  
}