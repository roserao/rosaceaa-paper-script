#!/usr/bin/env Rscript
library(rosace)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)

##### parse argument #####
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "Any of the valid data code. see data folder.")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data

##### key ##### 
key <- "1SM73"

ddir <- file.path("results_oct/models", data)
sdir <- file.path("results_oct/analysis", data)
if (!dir.exists(sdir)) {
  dir.create(sdir, recursive = TRUE)
}

sdir_table <- file.path(sdir, "table")
sdir_plot <- file.path(sdir, "plot")
if (!dir.exists(sdir_table)) {
  dir.create(sdir_table, recursive = TRUE)
}
if (!dir.exists(sdir_plot)) {
  dir.create(sdir_plot, recursive = TRUE)
}

########## 0 load full data ##########
load(file.path("results", "analysis", "oct1_sm73", "df_var_full.rda"))
df_var <- df_var3 
rm(df_var1, df_var2, df_var3)

########## 1 run time and cpu ##########
df_time <- rbind(
  read_tsv(file.path(ddir, "model1", "benchmark.txt")) %>% mutate(model = 1),
  read_tsv(file.path(ddir, "model2", "benchmark.txt")) %>% mutate(model = 2),
  read_tsv(file.path(ddir, "model3", "benchmark.txt")) %>% mutate(model = 3)
)
write_tsv(df_time, file = file.path(sdir_table, "benchmark.tsv"))

########## 2 combine object ##########
##### model 1 #####
rosace <- readRDS(file.path(ddir, "model1", "rosace_eval.rds"))
##### model 2 #####
rosace2 <- readRDS(file.path(ddir, "model2", "rosace_eval.rds"))
rosace <- AddScoreData(rosace, rosace2@scores[[1]])
##### model 3 #####
rosace3 <- readRDS(file.path(ddir, "model3", "rosace_eval.rds"))
rosace <- AddScoreData(rosace, rosace3@scores[[1]])
##### check #####
rm(rosace2, rosace3)
names(rosace@scores)

########## 3 check diagnosis ##########
df_diag <- rbind(
  data.frame(rosace@scores[[1]]@misc$diags) %>% mutate(chain = 1:4, model = 1),
  data.frame(rosace@scores[[2]]@misc$diags) %>% mutate(chain = 1:4, model = 2),
  data.frame(rosace@scores[[3]]@misc$diags) %>% mutate(chain = 1:4, model = 3)
)
write_tsv(df_diag, file = file.path(sdir_table, "mcmc_diags.tsv"))

########## 4 Output Score/Position ##########
##### model 1 #####
df_list1 <- OutputScore(rosace, name = paste(key, "ROSACE1", sep = "_"),
                        pos.info = TRUE,
                        sig.test = 0.05)
df_var1 <- df_list1$df_variant
df_pos1 <- df_list1$df_position
rm(df_list1)
##### model 2 #####
df_list2 <- OutputScore(rosace, name = paste(key, "ROSACE2", sep = "_"),
                        pos.info = TRUE, blosum.info = TRUE,
                        sig.test = 0.05)
df_var2 <- df_list2$df_variant
df_pos2 <- df_list2$df_position
rm(df_list2)
##### model 3 #####
df_list3 <- OutputScore(rosace, name = paste(key, "ROSACE3", sep = "_"),
                        pos.info = TRUE, blosum.info = TRUE, pos.act.info = TRUE,
                        sig.test = 0.05)
df_var3 <- df_list3$df_variant
df_pos3 <- df_list3$df_position
rm(df_list3)
##### save object #####
save(df_var1, df_var2, df_var3, file = file.path(sdir, "df_var_full.rda"))
save(df_pos1, df_pos2, df_pos3, file = file.path(sdir, "df_pos_full.rda"))

# load(file.path(sdir, "df_var_full.rda"))
# load(file.path(sdir, "df_pos_full.rda"))

########## 5 r-squared prediction value ##########
# df_nu_3
df_nu <- df_var3 %>% 
  select(blosum_score, stop, nu_mean, nu_sd) %>% unique() %>% 
  arrange(blosum_score, desc(stop)) 
write_tsv(df_nu, file.path(sdir_table, "df_nu3.tsv"))

# df_var_full
df_var_random <- df_var %>% 
  select(1:10, blosum_score) %>%
  filter(!ctrl, !stop) %>% 
  left_join(df_pos3 %>% select(pos, phi_mean, sigma2_mean, rho_mean)) %>%
  left_join(df_nu %>% select(blosum_score, stop, nu_mean)) %>%
  left_join(df_var3 %>% select(variants, mean_sample = mean))
df_var_random <- df_var_random %>% mutate(pred = phi_mean + rho_mean * nu_mean)
write_tsv(df_var_random , file.path(sdir_table, "df_var_full.tsv"))

# fit model r squared
fit1 <- summary(lm(mean ~ pred, data = df_var_random)) # 0.6511
fit2 <- summary(lm(mean ~ pred, data = df_var_random %>% filter(is.na(mean_sample))))
df_rsquared <- data.frame(
  r.squared_all = fit1$r.squared, 
  adj.r.squared_all = fit1$adj.r.squared,
  r.squared_missing = fit2$r.squared,
  adj.r.squared_missing = fit2$adj.r.squared
)
write_tsv(df_rsquared, file.path(sdir_table, "df_rsquared.tsv"))


######### finale ##########
saveRDS(rosace, file = file.path(sdir, "rosace_eval_full.rds"))
