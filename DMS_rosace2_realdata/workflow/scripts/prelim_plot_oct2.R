#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)

oct_list <- read_table("config/oct1_list.txt", col_names = "run_ID")


df_rsquared <- data.frame()
for (i in 1:nrow(oct_list)) {
    key <- oct_list$run_ID[i]
    ddir <- file.path("results_oct", "analysis", key)

    df_rsquared_sub <- read_tsv(file.path(ddir, "table", "df_rsquared.tsv"))
    df_rsquared <- rbind(df_rsquared, df_rsquared_sub %>% mutate(run_ID = key))
}
write_tsv(df_rsquared, file = "results_oct/summary/df_rsquared.tsv")

