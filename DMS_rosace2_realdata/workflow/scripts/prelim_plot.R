#!/usr/bin/env Rscript
library(rosace)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)

##### helper function #####
compute_rank_syn <- function(mean, lfsr, ctrl, rank.lfsr = TRUE, resolution = 100) {
  df <- data.frame(mean = mean, lfsr = lfsr, ctrl = ctrl)
  ### rank variants 
  if (rank.lfsr) {
    df <- df %>% arrange(lfsr, desc(abs(mean)))
  } else {
    df <- df %>% arrange(desc(abs(mean)), lfsr)
  }
  ### compute percentage of synonymous mutation called
  n_syn <- sum(df$ctrl)
  cutoff <- floor(nrow(df)/resolution * (1:resolution))
  result <- data.frame(rank = 1:resolution, fdr = -0.1)
  for (i in 1:resolution) {
    df_sig <- df[1:cutoff[i], ] %>% drop_na()
    result[i, 2] <- sum(df_sig$ctrl, na.rm = T)/n_syn
  }
  return(result)
}

##### parse argument #####
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "Any of the valid data code. see data folder."),
    make_option(c("-f", "--domain"), type = "logical", default = FALSE, 
        help = "whether to process domainase data"),
    make_option(c("-g", "--met"), type = "logical", default = FALSE, 
        help = "whether to process met2 data")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data
flag_domain <- opt$domain
flag_met <- opt$met

##### key ##### 
if (flag_domain) {
  key <- data
} else if (flag_met) {
  key <- substr(data, 6, nchar(data))
} else {
  if (data == "oct1_sm73") {
    key <- "1SM73"
  } else if (data == "slco1b1") {
    key <- "SLCO1B1.ctrl"
  } else if (data == "slco1b1_ros") {
    key <- "SLCO1B1.ros"
  } else if (data == "metex14") {
    key <- "MET1Ex14"
  } else if (data == "abcg2_dox") {
    key <- "ABCG2.dox"
  } else if (data == "abcg2_mtx") {
    key <- "ABCG2.mtx"
  } else if (data == "abcg2_sn38") {
    key <- "ABCG2.sn38"
  } else if (data == "msh2") {
    key <- "MSH2"
  } else if (data == "grb2") {
    key <- "GRB2"
  } else if (data == "grb2_gab2") {
    key <- "GRB2.gab2"
  } else if (data == "psd95") {
    key <- "PSD95"
  } else if (data == "psd95_cript") {
    key <- "PSD95.cript"
  } else if (data == "SPG1_Olson") {
    key <- "SPG1_Olson"
  } else if (data == "RNC") {
    key <- "RNC"
  } else if (data == "CD19") {
    key <- "CD19"
  } else if (data == "CAS9_pos") {
    key <- "CAS9_pos"
  } else if (data == "CAS9_neg") {
    key <- "CAS9_neg"
  } else if (data == "CAR11_lof") {
    key <- "CAR11_lof"
  } else if (data == "CAR11_gof") {
    key <- "CAR11_gof"
  } else {
    stop("data option not available.")
  }
}

########## 0 working directory ##########
if (flag_domain) {
  ddir <- file.path("results_domainase/models", data)
  sdir <- file.path("results_domainase/analysis", data)
} else if (flag_met) {
  ddir <- file.path("results_met/models", data)
  sdir <- file.path("results_met/analysis", data)
} else {
  ddir <- file.path("results/models", data)
  sdir <- file.path("results/analysis", data)
}
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

########## 5 beta value ##########
##### 5.1 beta correlation plot #####
df_beta <- data.frame(
  beta_model1 = df_var1$mean,
  beta_model2 = df_var2$mean,
  beta_model3 = df_var3$mean
)
p1 <- ggplot(df_beta, aes(beta_model1, beta_model2)) +
  geom_point(size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p2 <- ggplot(df_beta, aes(beta_model2, beta_model3)) +
  geom_point(size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p3 <- ggplot(df_beta, aes(beta_model3, beta_model1)) +
  geom_point(size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p <- plot_grid(p1, p2, p3, nrow = 1)
rm(p1, p2, p3)
# p 
save_plot(filename = file.path(sdir_plot, "beta_corr_models.png"), 
          plot = p, base_height = 5, base_width = 15)
##### 5.2 beta histogram #####
p <- ggplot(df_beta %>% pivot_longer(everything()), aes(value)) +
  geom_histogram(bins = 50, color = "grey") +
  facet_wrap(vars(name), ncol = 1) +
  theme_classic()
save_plot(filename = file.path(sdir_plot, "beta_histogram.png"), 
          plot = p, base_height = 7, base_width = 10)
rm(df_beta, p)
##### 5.3 rank variants (synonymous mutation) #####
df_rank_syn <- rbind(
  compute_rank_syn(mean = df_var1$mean, lfsr = df_var1$lfsr, ctrl = df_var1$ctrl,
                   rank.lfsr = TRUE, resolution = 100) %>% mutate(model = "model1"),
  compute_rank_syn(mean = df_var1$mean, lfsr = df_var2$lfsr, ctrl = df_var2$ctrl,
                   rank.lfsr = TRUE, resolution = 100) %>% mutate(model = "model2"),
  compute_rank_syn(mean = df_var1$mean, lfsr = df_var3$lfsr, ctrl = df_var3$ctrl,
                   rank.lfsr = TRUE, resolution = 100) %>% mutate(model = "model3")
)
p <- ggplot(df_rank_syn, aes(x = rank/100, y = fdr)) +
  geom_step(aes(color = model, linetype = model), 
            size = 0.5, alpha = 0.9) +
  labs(x = "ranked variants by hypothesis testing",
       y = "% of synonymous called",
       color = NULL, linetype = NULL)
save_plot(filename = file.path(sdir_plot, "rank_syn_fdr.png"), 
          plot = p, base_height = 4, base_width = 6)
rm(df_rank_syn, p)

########## 6 phi value ##########
##### 6.0 empirical phi #####
empr_phi1 <- (df_var1 %>% group_by(pos, ctrl, stop) %>% 
                summarise(empr_phi1 = mean(mean)) %>%
                filter(ctrl == FALSE, stop == FALSE))[["empr_phi1"]]
empr_phi2 <- (df_var2 %>% group_by(pos, ctrl, stop) %>% 
                summarise(empr_phi2 = mean(mean)) %>%
                filter(ctrl == FALSE, stop == FALSE))[["empr_phi2"]]
empr_phi3 <- (df_var3 %>% group_by(pos, ctrl, stop) %>% 
                summarise(empr_phi3 = mean(mean)) %>%
                filter(ctrl == FALSE, stop == FALSE))[["empr_phi3"]]
df_phi <- data.frame(
  phi_model1 = df_pos1$phi_mean,
  phi_model2 = df_pos2$phi_mean,
  phi_model3 = df_pos3$phi_mean,
  phi_empr1 = empr_phi1,
  phi_empr2 = empr_phi2,
  phi_empr3 = empr_phi3
)
rm(empr_phi1, empr_phi2, empr_phi3)
##### 6.1 phi model plot #####
p1 <- ggplot(df_phi, aes(phi_model1, phi_model2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p2 <- ggplot(df_phi, aes(phi_model2, phi_model3)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p3 <- ggplot(df_phi, aes(phi_model3, phi_model1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p <- plot_grid(p1, p2, p3, nrow = 1)
rm(p1, p2, p3)
save_plot(filename = file.path(sdir_plot, "phi_corr_models.png"), 
          plot = p, base_height = 5, base_width = 15)
##### 6.2 phi shrinkage plot #####
p1 <- ggplot(df_phi, aes(phi_empr1, phi_model1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p2 <- ggplot(df_phi, aes(phi_empr2, phi_model2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p3 <- ggplot(df_phi, aes(phi_empr3, phi_model3)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p <- plot_grid(p1, p2, p3, nrow = 1)
rm(p1, p2, p3)
save_plot(filename = file.path(sdir_plot, "phi_corr_empr.png"), 
          plot = p, base_height = 5, base_width = 15)
rm(df_phi, p)

########## 7 position-level estimate ##########
##### 7.1 phi and sigma2 #####
p1 <- ggplot(df_pos1, aes(phi_mean, sigma2_mean)) +
  geom_point()
p2 <- ggplot(df_pos2, aes(phi_mean, sigma2_mean)) +
  geom_point()
p3 <- ggplot(df_pos3, aes(phi_mean, sigma2_mean)) +
  geom_point()
p <- plot_grid(p1, p2, p3, nrow = 1)
rm(p1, p2, p3)
save_plot(filename = file.path(sdir_plot, "pos_phi_sigma2_scatter.png"), 
          plot = p, base_height = 5, base_width = 15)
##### 7.2 model 3: phi, rho, sigma 3 way #####
p1 <- ggplot(df_pos3, aes(phi_mean, rho_mean)) +
  geom_point()
p2 <- ggplot(df_pos3, aes(phi_mean, sigma2_mean)) +
  geom_point()
p3 <- ggplot(df_pos3, aes(rho_mean, sigma2_mean)) +
  geom_point()
p <- plot_grid(p1, p2, p3, nrow = 1)
rm(p1, p2, p3)
save_plot(filename = file.path(sdir_plot, "pos_model3_scatter.png"), 
          plot = p, base_height = 5, base_width = 15)

########## 8 nu value estimate ##########
########## nu value ##########
##### 8.1 summarise nu into a table #####
if (sum(df_var3$stop) == 0){
  df_nu <- df_var3 %>% 
    select(blosum_score, stop, nu_mean, nu_sd) %>% unique() %>% 
    arrange(blosum_score) %>%
    mutate(n = as.numeric(table(df_var3$blosum_score)))
} else {
  if (table(df_var3$blosum_score, df_var3$stop)[1, 1] == 0) {
    df_nu <- df_var3 %>% 
      select(blosum_score, stop, nu_mean, nu_sd) %>% unique() %>% 
      arrange(blosum_score, desc(stop)) %>%
      mutate(n = c(sum(df_var3$stop), as.numeric(table(df_var3$blosum_score, df_var3$stop)[-1, 1])))
  } else {
    df_nu <- df_var3 %>% 
      select(blosum_score, stop, nu_mean, nu_sd) %>% unique() %>% 
      arrange(blosum_score, desc(stop)) %>%
      mutate(n = c(sum(df_var3$stop), as.numeric(table(df_var3$blosum_score, df_var3$stop)[, 1])))
  }
}
write_tsv(df_nu, file = file.path(sdir_table, "df_nu3.tsv"))
if (sum(df_var2$stop) == 0){
  df_nu2 <- df_var2 %>% 
    select(blosum_score, stop, nu_mean, nu_sd) %>% unique() %>% 
    arrange(blosum_score) %>%
    mutate(n = as.numeric(table(df_var2$blosum_score)))
} else {
  if (table(df_var2$blosum_score, df_var2$stop)[1, 1] == 0) {
    df_nu2 <- df_var2 %>% 
      select(blosum_score, stop, nu_mean, nu_sd) %>% unique() %>% 
      arrange(blosum_score, desc(stop)) %>%
      mutate(n = c(sum(df_var2$stop), as.numeric(table(df_var2$blosum_score, df_var2$stop)[-1, 1])))
  } else {
    df_nu2 <- df_var2 %>% 
      select(blosum_score, stop, nu_mean, nu_sd) %>% unique() %>% 
      arrange(blosum_score, desc(stop)) %>%
      mutate(n = c(sum(df_var2$stop), as.numeric(table(df_var2$blosum_score, df_var2$stop)[, 1])))
  }
}
write_tsv(df_nu2, file = file.path(sdir_table, "df_nu2.tsv"))
##### 8.2 plot nu vlaue #####
p1 <- ggplot(df_var1 %>% mutate(residual = mean - phi_mean),
             aes(y = residual, x = as.factor(df_var3$blosum_score))) +
  geom_boxplot() +
  geom_jitter(size = 0.01, alpha = 0.5, color = "grey") +
  labs(x = "blosum", y = "beta - phi")
p2 <- ggplot(df_nu, aes(x = as.factor(blosum_score))) +
  geom_bar(aes(y = nu_mean), stat = "identity", 
           fill = "skyblue", alpha = 0.7) +
  geom_point(aes(y = nu_mean), colour="orange", alpha = 0.9) +
  geom_errorbar(aes(ymin = nu_mean - nu_sd, ymax = nu_mean + nu_sd), 
                width=0.3, colour="orange", alpha=0.9, size=1) +
  labs(x = "blosum", y = "nu", title = "model3") +
  geom_text(aes(y = max(nu_mean), label = n), size = 2, vjust = -1)
p3 <- ggplot(df_nu2, aes(x = as.factor(blosum_score))) +
  geom_bar(aes(y = nu_mean), stat = "identity", 
           fill = "skyblue", alpha = 0.7) +
  geom_point(aes(y = nu_mean), colour="orange", alpha = 0.9) +
  geom_errorbar(aes(ymin = nu_mean - nu_sd, ymax = nu_mean + nu_sd), 
                width=0.3, colour="orange", alpha=0.9, size=1) +
  labs(x = "blosum", y = "nu", title = "model2") +
  geom_text(aes(y = max(nu_mean), label = n), size = 2, vjust = -1)
p <- plot_grid(p1, p2, p3, nrow = 3)
save_plot(filename = file.path(sdir_plot, "nu_barplot.png"), 
          plot = p, base_height = 14, base_width = 8)
rm(p1, p2, p3, p)

########## 9 heatmap and some realistic thing ##########
npos <- 150
### 9.1 Heatmap: model 3 
scoreHeatmap(data = df_var3,
             type.col = "ctrl",
             ctrl.name = TRUE, # the control mutation name
             pos.col = "pos",
             wt.col = "wt",
             mut.col = "mut",
             score.col = "mean",
             npos = npos,
             savedir = sdir_plot, 
             name = "heatmap_model3",
             savepdf = FALSE,
             show = FALSE)
### 9.2 Heatmap: model 3
scoreVlnplot(data = df_var3,
             pos.col = "pos",
             wt.col = "wt",
             score.col = "mean",
             jitter = TRUE,
             npos = npos,
             pos.step = 1,
             ht = 10,
             wd = 20,
             ncol = 1,
             savedir = sdir_plot,
             name = "violin_dist_model3", 
             savepdf = FALSE, 
             show = FALSE)
### 9.3 Phi, Rho, line plot
df_pos3_plot <- df_var3 %>% filter(!ctrl) %>% 
  select(pos, wt, phi_mean, sigma2_mean, rho_mean, rho_sd) %>% unique()
scoreLnplot(data = df_pos3_plot,
            savedir = sdir_plot,
            pos.col = "pos",
            wt.col = "wt",
            mean.col = "phi_mean",
            sd.col = "sigma2_mean",
            npos = npos,
            ht = 20,
            wd = 20,
            name = "line_phi_model3",
            savepdf = FALSE, 
            show = FALSE)
scoreLnplot(data = df_pos3_plot,
            savedir = sdir_plot,
            pos.col = "pos",
            wt.col = "wt",
            mean.col = "rho_mean",
            sd.col = "rho_sd",
            npos = npos,
            ht = 20,
            wd = 20,
            name = "line_rho_model3",
            savepdf = FALSE, 
            show = FALSE)

######### finale ##########
saveRDS(rosace, file = file.path(sdir, "rosace_eval_full.rds"))
