#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)

domain_list <- read_table("config/domain_list.txt", col_names = "dom_ID")
domain_list <- domain_list %>% 
    separate(dom_ID, c("gene", "family", "fpos"), remove = FALSE, sep = "_")  


########## nu value PCA (cluster domain) ##########

df_nu3 <- data.frame()
for (i in nrow(domain_list)) {
    key <- domain_list$dom_ID[i]
    ddir <- file.path("results_domainase", "analysis", key)

    df_nu3_sub <- read_tsv(file.path(ddir, "table", "df_nu3.tsv"))
    df_nu3 <- rbind(df_nu3, df_nu3_sub %>% mutate(dom_ID = key))
}

df_nu3 <- df_nu3 %>%
    mutate(lfsr.neg = stats::pnorm(0, mean = nu_mean, sd = nu_sd, lower.tail = FALSE),
           lfsr.pos = stats::pnorm(0, mean = nu_mean, sd = nu_sd, lower.tail = TRUE),
           lfsr = min(.data$lfsr.neg, .data$lfsr.pos))

df_nu3_wide <- df_nu3 %>% 
    filter(!stop) %>%
    select(blosum_score, nu_mean, dom_ID) %>%
    pivot_wider(id_cols = "dom_ID", names_from = "blosum_score", values_from = "nu_mean") 

