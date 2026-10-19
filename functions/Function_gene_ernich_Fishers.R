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
  
  source("/Users/ana207/Desktop/R_Projects/catlas_enrichment/Function_Pairwise_fishers.R")


## Arguments:

#      1. wide_df: A wide format df--usually output from gene_enrich_wide_df function
#                  Should have a column named "marker_gene" a binary column that says if it is a marker or not

## function:

  enrichment_gene_fishers <- function(wide_df= df_1){
    
    ## looping through each variable name 
    var_names <- setdiff(colnames(wide_df), c("Gene", "marker_gene", "NA"))
    
    # performing 
    loop_results <- list()
    
    for (x in var_names) {
      loop_results[[x]] <- pairwise_fishers(df_data = wide_df,
                                            col_1 = "marker_gene",
                                            col_2 = x,
                                            col_name = x)
      
    }
    
    # binding rows and FDR
    fishers_res <- bind_rows(loop_results) %>%
      mutate(adj_p_val = p.adjust(p_val, method = "BH"))
    
    return(fishers_res)
    
  }