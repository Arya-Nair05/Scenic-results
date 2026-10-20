##----------------------------------------------------------------------------------
## Function to return a wider format data frame for doing multiple fishers enrichment test
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

## Arguments:
  
#      1. data: Data frame with all the marker peaks
#      2. enrich_df: Data frame of peaks/TFs/ with which we need to find enrichments
#                    Should have columns named chr, start, end
#      3. col_name: Name of the column in enrich_df that has the variable name, egs: cell_type/tf_name etc
#      4. padj_thresh: Adjusted p value threshold for selecting marker peaks
#      5. avg_log2_thresh: Avg log2FC threshold for selecting markers
#      6. category: Specify if we need enrichments for upregulated markers/downregulated markers
  
## function to perform multiple enrichments:
enrichment_wide_df <- function(data_df = df_1,
                               enrich_df= df_2,
                               col_name= NULL,
                               padj_thresh = 1e-5,
                               avg_log2_thresh= 0.5,
                               category= "Upregulated"){
  

  ## making granges of this df and splitting by cell type:
    if(is.null(col_name) == TRUE){
      stop("Column name of enrich_df having the variable names not specified")
    } 
  
  enrich_gr <- GRanges(seqnames = enrich_df$chr,
                    ranges = IRanges(start = enrich_df$start, end = enrich_df$end),
                     var_name= enrich_df[[col_name]])
  
  print(length(enrich_gr))
  
  enrich_gr_split <- split(enrich_gr, f = mcols(enrich_gr)$var_name)

# making a subset of the marker peaks:
  colnames(data_df)[1] <- "peak"  
  
  if(category =="Upregulated"){
    df_subset <- data_df %>%
      filter(p_val_adj < padj_thresh & avg_log2FC > avg_log2_thresh)
    
  } else {
    df_subset <- data_df %>%
      filter(p_val_adj < padj_thresh & avg_log2FC < avg_log2_thresh)
  }
    

## finding overlap of the marker peaks with each of the cell type CAtlas list:
  data_df <- data_df %>%
      separate(peak, into = c("chrom", "start", "end"), sep = "-", remove = FALSE)
    
  df_GR <- makeGRangesFromDataFrame(data_df, 
                                    seqnames.field = "chrom", 
                                    start.field = "start",
                                    end.field = "end")

  overlap_list <- lapply(enrich_gr_split, function(x){
    
    df <- findOverlapPairs(query = df_GR, subject = x)
    first_gr <- S4Vectors::first(df)
    second_gr <- S4Vectors::second(df)
    
    res_df <- data.frame(first=
                         paste(seqnames(first_gr),
                               start(first_gr),
                               end(first_gr),
                               sep = "-"),
                       second=
                         paste(seqnames(second_gr),
                               start(second_gr),
                               end(second_gr),
                               sep="-"))
    
    colnames(res_df) <- c("peak", "overlap_peak")
    return(res_df)
  })
  overlap_df <- bind_rows(overlap_list, .id = "var_name")

## making a master df
  enrich_df <- data_df %>%
    select(1) %>%
    mutate(marker_peak = peak %in% df_subset$peak)
  
  enrich_df <- left_join(enrich_df,overlap_df, by = "peak")

## making wider format df
  wide_df <- enrich_df %>%
    mutate(marker_peak = marker_peak*1) %>%
    mutate(overlap_peak = ifelse(!is.na(overlap_peak),1,0)) %>%
    distinct() %>% # if a peak overlapped with multiple peaks in one cell type then they need to be collapsed so that it can be added to the wider df
    pivot_wider(names_from = var_name,
                values_from = overlap_peak,
                values_fill = 0)

  return(wide_df)
  
}