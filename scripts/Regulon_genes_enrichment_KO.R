##-------------------------------------------------------------------------------------------------
## Enrichment of SC/DE/S2 marker genes with the target genes of each regulon
##-------------------------------------------------------------------------------------------------

## libraries
  library(tidyverse)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(GenomicRanges)
  library(ggthemes)
  
  set.seed(1234)
  
  source("/Users/ana207/Desktop/R_Projects/Scenic-results/functions/marker_gene_enrich_wide_df.R")
  source("/Users/ana207/Desktop/R_Projects/Scenic-results/functions/markers_enrich_fishers.R")
  source("/Users/ana207/Desktop/R_Projects/Scenic-results/functions/pairwise_fishers.R")
  
##-------------------------------------------------------------------------------------------------
## SC stage:
##-------------------------------------------------------------------------------------------------
## loading the regulons
  sc_reg <- read.table("data/D0_eRegulon_direct.tsv", sep = "\t", header = TRUE)

## loading the Bulk DEGs table:
  degs <- read.delim("data/bulk_rna_degs.tsv", sep = "\t", header = TRUE)
  
## filtering the DEGs per stage:
  degs_list <- list()
  for (stage in c("D0","DE","S2")){
    
    df_name <- paste0("deg_",stage)
    degs_list[[df_name]] <- degs %>%
      filter(Stage== stage) %>%
      select(Stage, GeneName,dsig)
  }
  list2env(degs_list, envir = .GlobalEnv)
  
## filtering the DEGs upregulated per stage:
  degs_list <- list()
  for (stage in c("D0","DE","S2")){
    
    df_name <- paste0("up_degs_",stage)
    degs_list[[df_name]] <- degs %>%
      filter(Stage== stage & dsig==1) %>%
      select(Stage, GeneName)
  }
  list2env(degs_list, envir = .GlobalEnv)
  
## SC stage:
  sc_wide_df <- gene_enrich_wide_df(data_df = deg_D0,
                                    enrich_df = sc_reg,
                                    col_name = "eRegulon_name",
                                    max_triplet_rank = 20000,
                                    category = "Up")

## Fishers
  sc_gene_fisher <- enrichment_fishers(wide_df = sc_wide_df, 
                                       category = "Gene")
  sc_gene_fisher$stage <- "SC"
  
## Plotting:
   sc_gene_fisher <- sc_gene_fisher[order(sc_gene_fisher$adj_p_val,decreasing = TRUE),]
   sc_gene_fisher$cell_type <- factor(sc_gene_fisher$cell_type, levels = unique(sc_gene_fisher$cell_type))
   
   ggplot(sc_gene_fisher, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
     geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
     theme_few(base_size = 13)+
     scale_fill_viridis()+
     geom_text(aes(label = round(odds_ratio,1)),size = 3)+
     theme(axis.text.x = element_text(size = 7))+
     labs(title = "SC:Enrichment of target genes of regulons with KO marker genes",
          x= "",
          y="Regulons")
 
##-------------------------------------------------------------------------------------------------
## DE stage:
##-------------------------------------------------------------------------------------------------
## loading the regulons
  de_reg <- read.table("data/DE_eRegulon_direct.tsv", sep = "\t", header = TRUE)
   
## wide format:
   de_wide_df <- gene_enrich_wide_df(data_df = deg_DE,
                                     enrich_df = de_reg,
                                     col_name = "eRegulon_name",
                                     max_triplet_rank = 20000,
                                     category = "Up")
   
   ## Fishers
   de_gene_fisher <- enrichment_fishers(wide_df = de_wide_df,
                                        category = "Gene")
   de_gene_fisher$stage <- "DE"
   
   ## Plotting:
   de_gene_fisher <- de_gene_fisher[order(de_gene_fisher$adj_p_val,decreasing = TRUE),]
   de_gene_fisher$cell_type <- factor(de_gene_fisher$cell_type, levels = unique(de_gene_fisher$cell_type))
   
   ggplot(de_gene_fisher, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
     geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
     theme_few(base_size = 13)+
     scale_fill_viridis()+
     geom_text(aes(label = round(odds_ratio,1)),size = 3)+
     theme(axis.text.x = element_text(size = 7))+
     labs(title = "DE:Enrichment of target genes of regulons with KO marker genes",
          x= "",
          y="Regulons")
   
##-------------------------------------------------------------------------------------------------
## DE stage:
##-------------------------------------------------------------------------------------------------
## loading the regulons
  s2_reg <- read.table("data/S2_eRegulon_direct.tsv", sep = "\t", header = TRUE)
   
## wide format:
  s2_wide_df <- gene_enrich_wide_df(data_df = deg_S2,
                                    enrich_df = s2_reg,
                                    col_name = "eRegulon_name",
                                    max_triplet_rank = 20000,
                                    category = "Up")
   
## Fishers
  s2_gene_fisher <- enrichment_fishers(wide_df = s2_wide_df,
                                       category = "Gene")
  s2_gene_fisher$stage <- "S2"
   
## Plotting:
  s2_gene_fisher <- s2_gene_fisher[order(s2_gene_fisher$adj_p_val,decreasing = TRUE),] 
  s2_gene_fisher_sub <- s2_gene_fisher %>%
    filter(odds_ratio>0)
  s2_gene_fisher_sub$cell_type <- factor(s2_gene_fisher_sub$cell_type, levels = unique(s2_gene_fisher_sub$cell_type))
   
  ggplot(s2_gene_fisher_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
     geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
     theme_few(base_size = 12)+
     scale_fill_viridis()+
     geom_text(aes(label = round(odds_ratio,1)),size = 3)+
     theme(axis.text.x = element_text(size = 5))+
     labs(title =  "S2:Enrichment of target genes of regulons with KO marker genes",
          x= "",
          y="Regulons")
  
  
##-------------------------------------------------------------------------------------------------
## Combining all stages
##-------------------------------------------------------------------------------------------------
  all_fishers <- rbind(sc_gene_fisher, de_gene_fisher, s2_gene_fisher)
  
  ## cleaning the df:
  all_fishers_sub <- all_fishers %>%
    filter(odds_ratio >0 & adj_p_val < 1e-5) %>%
    mutate(stage= factor(stage, levels=unique(stage))) %>%
    arrange(stage)
  
  all_fishers_sub$cell_type <- factor(all_fishers_sub$cell_type, levels = unique(all_fishers_sub$cell_type))
  
  ggplot(all_fishers_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity", width=0.5, color="black",linewidth = 0.01,height=1.5)+
    #geom_point(aes(size = -log10(adj_p_val), fill = odds_ratio),shape = 21) +
    #scale_size_continuous(name = "-log10(p_value)", range = c(1,10)) +
    theme_few(base_size = 12)+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)), size = 3)+
    theme(axis.text.y = element_text(size = 8))+
    labs(title = "Enrichment of Regulon associated genes with KO marker genes",
         x= "Stages",
         y="Regulons")

##------------------------------------------------------------------------------
## WT gained genes
##------------------------------------------------------------------------------
## SC stage:
  sc_wide_df_wt <- gene_enrich_wide_df(data_df = deg_D0,
                                    enrich_df = sc_reg,
                                    col_name = "eRegulon_name",
                                    max_triplet_rank = 20000,
                                    category = "Down")
  
  ## Fishers
  sc_gene_fisher_wt <- enrichment_gene_fishers(wide_df = sc_wide_df_wt)
  sc_gene_fisher_wt$stage <- "SC"
  
  ## Plotting:
  sc_gene_fisher_wt <- sc_gene_fisher_wt[order(sc_gene_fisher_wt$adj_p_val,decreasing = TRUE),]
  sc_gene_fisher_wt$cell_type <- factor(sc_gene_fisher_wt$cell_type, levels = unique(sc_gene_fisher_wt$cell_type))
  
  ggplot(sc_gene_fisher_wt, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
    theme_few(base_size = 13)+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)),size = 3)+
    theme(axis.text.x = element_text(size = 7))+
    labs(title = "SC:Enrichment of target genes of regulons with WT marker genes",
         x= "",
         y="Regulons")
  
##-------------------------------------------------------------------------------------------------
## DE stage:
##-------------------------------------------------------------------------------------------------
## loading the regulons
  de_reg <- read.table("data/DE_eRegulon_direct.tsv", sep = "\t", header = TRUE)
  
  ## wide format:
  de_wide_df_wt <- gene_enrich_wide_df(data_df = deg_DE,
                                    enrich_df = de_reg,
                                    col_name = "eRegulon_name",
                                    max_triplet_rank = 20000,
                                    category = "Down")
  
  ## Fishers
  de_gene_fisher_wt <- enrichment_gene_fishers(wide_df = de_wide_df_wt)
  de_gene_fisher_wt$stage <- "DE"
  
  ## Plotting:
  de_gene_fisher_wt <- de_gene_fisher_wt[order(de_gene_fisher_wt$adj_p_val,decreasing = TRUE),]
  de_gene_fisher_wt$cell_type <- factor(de_gene_fisher_wt$cell_type, levels = unique(de_gene_fisher_wt$cell_type))
  
  ggplot(de_gene_fisher_wt, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
    theme_few(base_size = 13)+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)),size = 3)+
    theme(axis.text.x = element_text(size = 7))+
    labs(title = "DE:Enrichment of target genes of regulons with WT marker genes",
         x= "",
         y="Regulons")
  
##-------------------------------------------------------------------------------------------------
## S2 stage:
##-------------------------------------------------------------------------------------------------
## loading the regulons
  s2_reg <- read.table("data/S2_eRegulon_direct.tsv", sep = "\t", header = TRUE)
  
  ## wide format:
  s2_wide_df_wt <- gene_enrich_wide_df(data_df = deg_S2,
                                    enrich_df = s2_reg,
                                    col_name = "eRegulon_name",
                                    max_triplet_rank = 20000,
                                    category = "Down")
  
  ## Fishers
  s2_gene_fisher_wt <- enrichment_gene_fishers(wide_df = s2_wide_df_wt)
  s2_gene_fisher_wt$stage <- "S2"
  
  ## Plotting:
  s2_gene_fisher_wt <- s2_gene_fisher_wt[order(s2_gene_fisher_wt$adj_p_val,decreasing = TRUE),] 
  
  s2_gene_fisher_wt_sub <- s2_gene_fisher_wt %>%
    filter(odds_ratio>0)
  s2_gene_fisher_wt_sub$cell_type <- factor(s2_gene_fisher_wt_sub$cell_type, levels = unique(s2_gene_fisher_wt_sub$cell_type))
  
  ggplot(s2_gene_fisher_wt_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
    theme_few(base_size = 12)+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)),size = 3)+
    theme(axis.text.x = element_text(size = 5))+
    labs(title =  "S2:Enrichment of target genes of regulons with WT marker genes",
         x= "",
         y="Regulons")
  
  
##-------------------------------------------------------------------------------------------------
## Combining all stages
##-------------------------------------------------------------------------------------------------
  all_fishers_wt <- rbind(sc_gene_fisher_wt, de_gene_fisher_wt, s2_gene_fisher_wt)
  
  ## cleaning the df:
  all_fishers_wt_sub <- all_fishers_wt %>%
    filter(odds_ratio >0 & adj_p_val < 1e-5) %>%
    mutate(stage= factor(stage, levels=unique(stage))) %>%
    arrange(stage)
  
  all_fishers_wt_sub$cell_type <- factor(all_fishers_wt_sub$cell_type, levels = unique(all_fishers_wt_sub$cell_type))
  
  ggplot(all_fishers_wt_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity", width=0.5, color="black",linewidth = 0.01,height=1.5)+
    #geom_point(aes(size = -log10(adj_p_val), fill = odds_ratio),shape = 21) +
    #scale_size_continuous(name = "-log10(p_value)", range = c(1,10)) +
    theme_few(base_size = 12)+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)), size = 3)+
    theme(axis.text.y = element_text(size = 8))+
    labs(title = "Enrichment of Regulon associated genes with WT marker genes",
         x= "Stages",
         y="Regulons")
  
##-------------------------------------------------------------------------------------------------
## Combining KO and WT
##-------------------------------------------------------------------------------------------------
## laoding the fishers results
  all_fishers_wt$stage <- paste0(all_fishers_wt$stage, "_WT")
  all_fishers$stage <- paste0(all_fishers$stage, "_KO")
  
  ko_wt_fishers <- rbind(all_fishers, all_fishers_wt)

  ko_wt_fishers$stage <- factor(ko_wt_fishers$stage, levels = c("SC_KO","SC_WT","DE_KO","DE_WT","S2_KO","S2_WT"))
  
  ko_wt_fishers_sub <- ko_wt_fishers %>%
    group_by(stage) %>%
    arrange(adj_p_val, decreasing = FALSE) %>%
    slice_head(n=10) %>%
    ungroup() 
  
  reg_sub <- ko_wt_fishers %>%
    group_by(stage) %>%
    arrange(adj_p_val, decreasing = FALSE) %>%
    slice_head(n=18) %>%
    ungroup() %>%
    pull(cell_type)
    
  ko_wt_fishers_sub <- ko_wt_fishers[ko_wt_fishers$cell_type %in% reg_sub, ] 
  
  ggplot(ko_wt_fishers_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity", width=1, color="black",linewidth = 0.01,height=1.2)+
    theme_few(base_size = 13)+
    scale_fill_viridis(trans="sqrt")+
    geom_text(aes(label = round(odds_ratio,1)), size = 3)+
    theme(axis.text.y = element_text(size = 7.5))+
    labs(title = "Enrichment of Regulon associated genes with DEGs",
         x= "",
         y="Regulons")
  
 
  