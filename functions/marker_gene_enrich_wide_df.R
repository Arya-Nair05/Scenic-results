##----------------------------------------------------------------------------------
## Function to return a wider format data frame for doing multiple fishers enrichment test
## This for doing enrichments of marker genes and not marker peaks
##----------------------------------------------------------------------------------

## libraries
  library(tidyverse)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(GenomicRanges)
  
  set.seed(1234)

## Arguments:

#      1. data: Data frame of DEGs per stage (should include all genes and not just the up genes)
#               Should have a column named Gene
#      2. enrich_df: Data frame of peaks/TFs/genes with which we need to find enrichments
#                    Should have a column named Gene
#      3. col_name: Name of the column in enrich_df that has the variable name, egs: regulon_name/cell_type/tf_name etc
#      4. category: Specify if we need enrichments for upregulated markers/downregulated markers
#      5. max_triplet_rank: The max number of regulons to be selecetd based on the triplet 
  
## function to perform multiple enrichments:
gene_enrich_wide_df <- function(data_df = df_1,
                                 enrich_df= df_2,
                                 col_name= NULL,
                                 max_triplet_rank= 20000,
                                 category= "Up"){
    
  
  ## making granges of this df and splitting by cell type:
    if(is.null(col_name) == TRUE){
      stop("Column name of enrich_df having the variable names not specified")
    } 
  
  ## splitting the regulon df:
    enrich_df_filt <- enrich_df %>%
      filter(triplet_rank < max_triplet_rank) %>%
      filter(rho_TF2G >0) %>%
      select(Gene, col_name)
    
    enrich_df_filt[[col_name]] <- factor(enrich_df_filt[[col_name]], 
                                        levels = unique(enrich_df_filt[[col_name]]))
    
    enrich_gr_split <- split(enrich_df_filt, f=enrich_df_filt$eRegulon_name)

  ## making a subset of the marker peaks:
    colnames(data_df)[2] <- "Gene"  
  
  ## finding all the up genes:
    
    if (category == "Up"){
      degs_filt <- data_df[data_df$dsig==1,]
      
    } else {
      degs_filt <- data_df[data_df$dsig==-1,]
    }
      
    
  ## finding overlap of the marker genes with each of the regulon target genes:
    for (reg in names(enrich_gr_split)){
      
      df <- enrich_gr_split[[reg]]
      df$up_deg <- df$Gene %in% degs_filt$Gene
      enrich_gr_split[[reg]] <- df
      
    }
    
    overlap_df <- bind_rows(enrich_gr_split, .id = "eRegulon_name") %>%
      distinct()
  
  ## making a master df
    deg_all <- data_df %>%
      select(-c(Stage,dsig)) %>%
      mutate(marker_gene= Gene %in% degs_filt$Gene) %>%
      distinct()
  
   enrich_df <- left_join(deg_all,overlap_df, by = "Gene")
  
  ## making wider format df
    gene_wide_df <- enrich_df %>%
      mutate(marker_gene= marker_gene*1,
             up_deg = up_deg*1) %>%
      pivot_wider(names_from = eRegulon_name,
                  values_from = up_deg,
                  values_fill = 0)
    
  return(gene_wide_df)
  
}