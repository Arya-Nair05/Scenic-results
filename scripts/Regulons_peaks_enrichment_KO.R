##-------------------------------------------------------------------------------------------------
## Enrichment of SC/DE/S2 marker peaks with regions associated with each regulons
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
  
  source("../../Analysis/smg_functions/Function_Overlaps_Output_Bed.R")
  source("../catlas_enrichment/Function_Pairwise_fishers.R")
  source("../catlas_enrichment/Function_enrichment_wide_df.R")
  source("../catlas_enrichment/Function_enrichment_fishers.R")
  
##-------------------------------------------------------------------------------------------------
## SC stage:
##-------------------------------------------------------------------------------------------------
## loading the regulons
  sc_reg <- read.table("D0_eRegulon_direct.tsv", sep = "\t", header = TRUE)
  
## selcting the top 20k regulon triplets and cleaning the df
  reg_num <- 20000
  
  sc_reg_filt <- sc_reg %>%
    filter(triplet_rank < reg_num) %>%
    filter(rho_TF2G >0) %>%
    select(Region, eRegulon_name) %>%
    separate(Region, into = c("chrom", "start", "end"), sep = ":|-") %>%
    mutate(start= as.numeric(start),
           end= as.numeric(end))
  
## loading the markers:
  sc_markers <- read.csv("../catlas_enrichment/SC_KO_vs_WT_markers.csv", header = TRUE)
  
## Fishers:
  # making a wide format df:
  sc_wide <- enrichment_wide_df(data_df = sc_markers,
                                enrich_df = sc_reg_filt,
                                col_name = "eRegulon_name",
                                padj_thresh = 1e-5,
                                avg_log2_thresh = 0.7,
                                category = "Upregulated")
  
  sc_fishers <- enrichment_fishers(wide_df = sc_wide)
  sc_fishers$stage <- "SC"
  
## heatmap and dotplot:
  sc_fishers <- sc_fishers[order(sc_fishers$adj_p_val, decreasing = TRUE),]
  sc_fishers$cell_type <- factor(sc_fishers$cell_type, levels = unique(sc_fishers$cell_type))
  sc_fishers_sub <- sc_fishers[sc_fishers$odds_ratio >0,]
  
  
  ggplot(sc_fishers_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
    theme_few()+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)))+
    theme(axis.text.x = element_text(size = 12))+
    labs(title = "SC:Enrichment of Regulon associated peaks with KO marker peaks",
         x= "")
 
##-------------------------------------------------------------------------------------------------
## DE stage
##-------------------------------------------------------------------------------------------------
## loading the regulons
  de_reg <- read.table("DE_eRegulon_direct.tsv", sep = "\t", header = TRUE)
  
## selcting the top 20k regulon triplets and removing regulons with negative TF-GENE rho value
  reg_num <- 20000
  
  de_reg_filt <- de_reg %>%
    filter(triplet_rank < reg_num) %>%
    filter(rho_TF2G >0) %>%
    select(Region, eRegulon_name) %>%
    separate(Region, into = c("chrom", "start", "end"), sep = ":|-") %>%
    mutate(start= as.numeric(start),
           end= as.numeric(end)) 
  
  ## loading the markers:
  de_markers <- read.csv("../catlas_enrichment/DE_KO_vs_WT_markers.csv", header = TRUE)
  
  ## Fishers:
  # making a wide format df:
  de_wide <- enrichment_wide_df(data_df = de_markers,
                                enrich_df = de_reg_filt,
                                col_name = "eRegulon_name",
                                padj_thresh = 1e-5,
                                avg_log2_thresh = 0.7,
                                category = "Upregulated")
  
  de_fishers <- enrichment_fishers(wide_df = de_wide)
  de_fishers$stage <- "DE"
  
  ## heatmap and dotplot:
  de_fishers <- de_fishers[order(de_fishers$adj_p_val, decreasing = TRUE),]
  de_fishers$cell_type <- factor(de_fishers$cell_type, levels = unique(de_fishers$cell_type))
  de_fishers_sub <- de_fishers[de_fishers$odds_ratio >0,]
  
  ggplot(de_fishers_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
    theme_few()+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)))+
    theme(axis.text.x = element_text(size = 12))+
    labs(title = "DE:Enrichment of Regulon associated peaks with KO marker peaks",
         x= "")
##-------------------------------------------------------------------------------------------------
## S2 stage
##-------------------------------------------------------------------------------------------------
## loading the regulons
  s2_reg <- read.table("S2_eRegulon_direct.tsv", sep = "\t", header = TRUE)
  
  ## selcting the top 20k regulon triplets and cleaning the df
  reg_num <- 20000
  
  s2_reg_filt <- s2_reg %>%
    filter(triplet_rank < reg_num) %>%
    filter(rho_TF2G >0) %>%
    select(Region, eRegulon_name) %>%
    separate(Region, into = c("chrom", "start", "end"), sep = ":|-") %>%
    mutate(start= as.numeric(start),
           end= as.numeric(end)) 
  
  ## loading the markers:
  s2_markers <- read.csv("../catlas_enrichment/S2_KO_vs_WT_markers.csv", header = TRUE)
  
  ## Fishers:
  # making a wide format df:
  s2_wide <- enrichment_wide_df(data_df = s2_markers,
                                enrich_df = s2_reg_filt,
                                col_name = "eRegulon_name",
                                padj_thresh = 1e-5,
                                avg_log2_thresh = 0.7,
                                category = "Upregulated")
  
  s2_fishers <- enrichment_fishers(wide_df = s2_wide)
  s2_fishers$stage <- "S2"
  
  ## heatmap and dotplot:
  s2_fishers <- s2_fishers[order(s2_fishers$adj_p_val, decreasing = TRUE),]
  s2_fishers$cell_type <- factor(s2_fishers$cell_type, levels = unique(s2_fishers$cell_type))
  s2_fishers_sub<- s2_fishers[s2_fishers$odds_ratio >0,]
  
  ggplot(s2_fishers_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
    theme_few()+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)))+
    theme(axis.text.x = element_text(size = 12))+
    labs(title = "S2:Enrichment of Regulon associated peaks with KO marker peaks",
         x= "")
  
##-------------------------------------------------------------------------------------------------
## Combining all stages
##-------------------------------------------------------------------------------------------------
  all_fishers <- rbind(sc_fishers, de_fishers, s2_fishers)
  
## cleaning the df:
  all_fishers_sub <- all_fishers %>%
    filter(odds_ratio >0 & adj_p_val < 0.005) %>%
    mutate(stage= factor(stage, levels=unique(stage))) %>%
    arrange(stage)
  
  all_fishers_sub$cell_type <- factor(all_fishers_sub$cell_type, levels = unique(all_fishers_sub$cell_type))
  
  ggplot(all_fishers_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity", width=0.8, color="black",linewidth = 0.01)+
    #geom_point(aes(size = -log10(adj_p_val), fill = odds_ratio),shape = 21) +
    #scale_size_continuous(name = "-log10(p_value)", range = c(1,10)) +
    theme_few(base_size = 14)+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)), size = 3)+
    theme(axis.text.y = element_text(size = 8))+
    labs(title = "Enrichment of Regulon associated peaks with KO marker peaks",
         x= "Stages",
         y="Regulons")
     
  
  
##-----------------------------------------------------------------------------
## WT marker peaks
##----------------------------------------------------------------------------
## SC stage:
## loading the regulons
  sc_reg <- read.table("D0_eRegulon_direct.tsv", sep = "\t", header = TRUE)
  
## selcting the top 20k regulon triplets and cleaning the df
  reg_num <- 20000
  
  sc_reg_filt <- sc_reg %>%
    filter(triplet_rank < reg_num) %>%
    filter(rho_TF2G >0) %>%
    select(Region, eRegulon_name) %>%
    separate(Region, into = c("chrom", "start", "end"), sep = ":|-") %>%
    mutate(start= as.numeric(start),
           end= as.numeric(end))
  
## loading the markers:
  sc_markers <- read.csv("../catlas_enrichment/SC_KO_vs_WT_markers.csv", header = TRUE)
  
## Fishers:
# making a wide format df:
  sc_wide_wt <- enrichment_wide_df(data_df = sc_markers,
                                enrich_df = sc_reg_filt,
                                col_name = "eRegulon_name",
                                padj_thresh = 1e-5,
                                avg_log2_thresh = 0.7,
                                category = "Down")
  
  sc_fishers_wt <- enrichment_fishers(wide_df = sc_wide_wt)
  sc_fishers_wt$stage <- "SC"
  
  ## heatmap and dotplot:
  sc_fishers_wt <- sc_fishers_wt[order(sc_fishers_wt$adj_p_val, decreasing = TRUE),]
  sc_fishers_wt$cell_type <- factor(sc_fishers_wt$cell_type, levels = unique(sc_fishers_wt$cell_type))
  sc_fishers_wt_sub <- sc_fishers_wt[sc_fishers_wt$odds_ratio >0,]
  
  
  ggplot(sc_fishers_wt_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
    theme_few()+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)))+
    theme(axis.text.x = element_text(size = 12))+
    labs(title = "SC:Enrichment of Regulon associated peaks with WT marker peaks",
         x= "")
  
##-------------------------------------------------------------------------------------------------
## DE stage
##-------------------------------------------------------------------------------------------------
## loading the regulons
  de_reg <- read.table("DE_eRegulon_direct.tsv", sep = "\t", header = TRUE)
  
## selcting the top 20k regulon triplets and removing regulons with negative TF-GENE rho value
  reg_num <- 20000
  
  de_reg_filt <- de_reg %>%
    filter(triplet_rank < reg_num) %>%
    filter(rho_TF2G >0) %>%
    select(Region, eRegulon_name) %>%
    separate(Region, into = c("chrom", "start", "end"), sep = ":|-") %>%
    mutate(start= as.numeric(start),
           end= as.numeric(end)) 
  
## loading the markers:
  de_markers <- read.csv("../catlas_enrichment/DE_KO_vs_WT_markers.csv", header = TRUE)
  
## Fishers:
  # making a wide format df:
  de_wide_wt <- enrichment_wide_df(data_df = de_markers,
                                enrich_df = de_reg_filt,
                                col_name = "eRegulon_name",
                                padj_thresh = 1e-5,
                                avg_log2_thresh = 0.7,
                                category = "Down")
  
  de_fishers_wt <- enrichment_fishers(wide_df = de_wide_wt)
  de_fishers_wt$stage <- "DE"
  
  ## heatmap and dotplot:
  de_fishers_wt <- de_fishers_wt[order(de_fishers_wt$adj_p_val, decreasing = TRUE),]
  de_fishers_wt$cell_type <- factor(de_fishers_wt$cell_type, levels = unique(de_fishers_wt$cell_type))
  de_fishers_wt_sub <- de_fishers_wt[de_fishers_wt$odds_ratio >0,]
  
  ggplot(de_fishers_wt_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
    theme_few()+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)))+
    theme(axis.text.x = element_text(size = 12))+
    labs(title = "DE:Enrichment of Regulon associated peaks with WT marker peaks",
         x= "")
##-------------------------------------------------------------------------------------------------
## S2 stage
##-------------------------------------------------------------------------------------------------
## loading the regulons
  s2_reg <- read.table("S2_eRegulon_direct.tsv", sep = "\t", header = TRUE)
  
## selcting the top 20k regulon triplets and cleaning the df
  reg_num <- 20000
  
  s2_reg_filt <- s2_reg %>%
    filter(triplet_rank < reg_num) %>%
    filter(rho_TF2G >0) %>%
    select(Region, eRegulon_name) %>%
    separate(Region, into = c("chrom", "start", "end"), sep = ":|-") %>%
    mutate(start= as.numeric(start),
           end= as.numeric(end)) 
  
## loading the markers:
  s2_markers <- read.csv("../catlas_enrichment/S2_KO_vs_WT_markers.csv", header = TRUE)
  
## Fishers:
  # making a wide format df:
  s2_wide_wt <- enrichment_wide_df(data_df = s2_markers,
                                enrich_df = s2_reg_filt,
                                col_name = "eRegulon_name",
                                padj_thresh = 1e-5,
                                avg_log2_thresh = 0.7,
                                category = "Down")
  
  s2_fishers_wt <- enrichment_fishers(wide_df = s2_wide_wt)
  s2_fishers_wt$stage <- "S2"
  
  ## heatmap and dotplot:
  s2_fishers_wt <- s2_fishers_wt[order(s2_fishers_wt$adj_p_val, decreasing = TRUE),]
  s2_fishers_wt$cell_type <- factor(s2_fishers_wt$cell_type, levels = unique(s2_fishers_wt$cell_type))
  s2_fishers_wt_sub <- s2_fishers_wt[s2_fishers_wt$odds_ratio >0,]
  
  ggplot(s2_fishers_wt_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity",width=0.5,color="black",linewidth = 0.01)+
    theme_few()+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)))+
    theme(axis.text.x = element_text(size = 12))+
    labs(title = "S2:Enrichment of Regulon associated peaks with WT marker peaks",
         x= "")
  
##-------------------------------------------------------------------------------------------------
## Combining all stages
##-------------------------------------------------------------------------------------------------
  all_fishers_wt <- rbind(sc_fishers_wt, de_fishers_wt, s2_fishers_wt)
  
  ## cleaning the df:
  all_fishers_wt_sub <- all_fishers_wt %>%
    filter(odds_ratio >0 & adj_p_val < 0.005) %>%
    mutate(stage= factor(stage, levels=unique(stage))) %>%
    arrange(stage)
  
  all_fishers_wt_sub$cell_type <- factor(all_fishers_wt_sub$cell_type, levels = unique(all_fishers_wt_sub$cell_type))
  
  ggplot(all_fishers_wt_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity", width=0.5, color="black",linewidth = 0.01)+
    #geom_point(aes(size = -log10(adj_p_val), fill = odds_ratio),shape = 21) +
    #scale_size_continuous(name = "-log10(p_value)", range = c(1,10)) +
    theme_few(base_size = 12)+
    scale_fill_viridis()+
    geom_text(aes(label = round(odds_ratio,1)), size = 3)+
    theme(axis.text.y = element_text(size = 8))+
    labs(title = "Enrichment of Regulon associated peaks with WT marker peaks",
         x= "Stages",
         y="Regulons")
  
##-------------------------------------------------------------------------------------------------
## Combining WT and KO
##-------------------------------------------------------------------------------------------------
 ## adding genotype info along with the stage
  all_fishers$stage <- paste0(all_fishers$stage, "_KO")
  all_fishers_wt$stage <- paste0(all_fishers_wt$stage, "_WT")
  
  regs_ko_wt <- rbind(all_fishers, all_fishers_wt)
  regs_ko_wt$stage <- factor(regs_ko_wt$stage, levels = c("SC_KO","SC_WT","DE_KO","DE_WT","S2_KO","S2_WT"))
  
  reg_ko_wt_sub <- regs_ko_wt %>%
    group_by(stage) %>%
    arrange(adj_p_val, decreasing = FALSE) %>%
    slice_head(n=10) %>%
    ungroup()
  
  regulons_sub <- regs_ko_wt %>%
    group_by(stage) %>%
    arrange(adj_p_val, decreasing = FALSE) %>%
    slice_head(n=18) %>%
    ungroup() %>%
    pull(cell_type)
  
  reg_ko_wt_sub <- regs_ko_wt%>%
    filter(cell_type %in% regulons_sub) %>%
    arrange(stage, desc(-(adj_p_val)))
  
  reg_ko_wt_sub$cell_type <- factor(reg_ko_wt_sub$cell_type, levels = unique(reg_ko_wt_sub$cell_type))
  
  ggplot(reg_ko_wt_sub, aes(x=stage, y = cell_type, fill = -log10(adj_p_val)))+
    geom_tile(stat = "identity", width=1, color="black",linewidth = 0.01)+
    #geom_point(aes(size = -log10(adj_p_val), fill = odds_ratio),shape = 21) +
    #scale_size_continuous(name = "-log10(p_value)", range = c(1,10)) +
    theme_few(base_size = 12)+
    scale_fill_viridis(trans="sqrt")+
    geom_text(aes(label = round(odds_ratio,1)), size = 3)+
    theme(axis.text.y = element_text(size = 8))+
    labs(title = "Enrichment of Regulon associated peaks the marker peaks",
         x= "",
         y="Regulons")
  