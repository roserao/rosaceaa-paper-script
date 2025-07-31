#!/usr/bin/env Rscript
library('rosace')
library('dplyr')
library('tidyr')
library('readr')
library('ggplot2')
library('cowplot')
source("workflow/scripts/analyze_rosace_utils.R")

##### parse argument #####
library('optparse')
option_list <- list(
    make_option(c("-m", "--mode"), type = "numeric", default = NULL, 
        help = "simualation mode: 1, 2, 3, 4, 5"),
    make_option(c("-d", "--data"), type = "numeric", default = NULL, 
        help = "data code: 1(abcg2), 2(oct1+sm73), 3(slco1b1+mtx), 4(abcg2+ros)"),
    make_option(c("-r", "--rep"), type = "numeric", default = NULL,
        help = "number of replicates"),
    make_option(c("-t", "--round"), type = "numeric", default = NULL, 
        help = "number of rounds"),
    make_option(c("-s", "--sim"), type = "numeric", default = NULL, 
        help = "simulation index")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data_mode <- opt$data
sim_mode <- opt$mode
n_rep <- opt$rep
n_round <- opt$round
idx_sim <- opt$sim

# data_mode <- 1
# sim_mode <- 1
# n_rep <- 3
# n_round <- 3
# idx_sim <- 2

##### 0 directory #####
ddir <-  
    file.path('results/models', 
              paste('data', data_mode, sep = ""),
              paste("growth_rep", n_rep, "_rd", n_round, "_", sim_mode, sep = ""),
              paste("sim", idx_sim, sep = ""))
sdir <-
    file.path('results/analysis', 
              paste('data', data_mode, sep = ""),
              paste("growth_rep", n_rep, "_rd", n_round, "_", sim_mode, sep = ""),
              paste("sim", idx_sim, sep = ""))  
if (!dir.exists(sdir)) {
  dir.create(sdir, recursive = TRUE)
}

########## read time and cpu usage ##########
df_time <- rbind(
  read_tsv(file.path(ddir, "model1", "benchmark.txt")) %>% mutate(model = 1),
  read_tsv(file.path(ddir, "model2", "benchmark.txt")) %>% mutate(model = 2),
  read_tsv(file.path(ddir, "model3", "benchmark.txt")) %>% mutate(model = 3)
)
write_tsv(df_time, file = file.path(sdir, "benchmark.tsv"))

########## 1 combine object ##########
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
# saveRDS(rosace, file = file.path(ddir, "rosace_eval_full.rds"))
# rosace <- readRDS(file.path(ddir, "rosace_eval_full.rds"))

########## 2 check diagnosis ##########
df_diag <- rbind(
  data.frame(rosace@scores[[1]]@misc$diags) %>% mutate(model = 1),
  data.frame(rosace@scores[[2]]@misc$diags) %>% mutate(model = 2),
  data.frame(rosace@scores[[3]]@misc$diags) %>% mutate(model = 3)
)
write_tsv(df_diag, file = file.path(sdir, "mcmc_diags.tsv"))

########## 3 Output Score/Position ##########
##### model 1 #####
df_list1 <- OutputScore(rosace, name = "simulation_ROSACE1", 
                        pos.info = TRUE,
                        sig.test = 0.05)
df_var1 <- df_list1$df_variant
df_pos1 <- df_list1$df_position
rm(df_list1)
##### model 2 #####
df_list2 <- OutputScore(rosace, name = "simulation_ROSACE2", 
                        pos.info = TRUE, blosum.info = TRUE,
                        sig.test = 0.05)
df_var2 <- df_list2$df_variant
df_pos2 <- df_list2$df_position
rm(df_list2)
##### model 3 #####
df_list3 <- OutputScore(rosace, name = "simulation_ROSACE3", 
                        pos.info = TRUE, blosum.info = TRUE, pos.act.info = TRUE,
                        sig.test = 0.05)
df_var3 <- df_list3$df_variant
df_pos3 <- df_list3$df_position
rm(df_list3)
##### save object #####
save(df_var1, df_var2, df_var3, file = file.path(sdir, "df_var_full.rda"))
save(df_pos1, df_pos2, df_pos3, file = file.path(sdir, "df_pos_full.rda"))

########## 4 beta value and lfsr (fdr and sensitivity) ##########
##### 4.1 beta correlation #####
df_corr <- 
  data.frame(
    model = 1:3,
    pearson = c(
      cor(df_var1$expected_effect, df_var1$mean, method = "pearson"),
      cor(df_var2$expected_effect, df_var2$mean, method = "pearson"),
      cor(df_var3$expected_effect, df_var3$mean, method = "pearson")
    ),
    spearman = c(
      cor(df_var1$expected_effect, df_var1$mean, method = "spearman"),
      cor(df_var2$expected_effect, df_var2$mean, method = "spearman"),
      cor(df_var3$expected_effect, df_var3$mean, method = "spearman")
    )
  )
write_tsv(df_corr, file = file.path(sdir, "beta_corr.tsv"))
rm(df_corr)
##### 4.2 beta correlation plot #####
df_beta <- data.frame(
  beta_true = df_var1$expected_effect,
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
save_plot(filename = file.path(sdir, "beta_corr_models.png"), 
          plot = p, base_height = 5, base_width = 15)
p1 <- ggplot(df_beta, aes(beta_true, beta_model1)) +
  geom_point(size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p2 <- ggplot(df_beta, aes(beta_true, beta_model2)) +
  geom_point(size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p3 <- ggplot(df_beta, aes(beta_true, beta_model3)) +
  geom_point(size = 0.5) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p <- plot_grid(p1, p2, p3, nrow = 1)
rm(p1, p2, p3)
# p
save_plot(filename = file.path(sdir, "beta_corr_expected.png"), 
          plot = p, base_height = 5, base_width = 15)
### histogram 
p <- ggplot(df_beta %>% pivot_longer(everything()), aes(value)) +
  geom_histogram(bins = 50, color = "grey") +
  facet_wrap(vars(name)) +
  theme_classic()
save_plot(filename = file.path(sdir, "beta_histogram.png"), 
          plot = p, base_height = 7, base_width = 10)
rm(df_beta, p)
##### 4.3 fdr and sensitivity #####
# FDR = FP/(FP+TP)
# Sensitivity = TP/(TP+FN)
expected_label <- df_var1$expected_label
expected_label[expected_label == "ctrl"] <- "Neutral"
df_rate <- data.frame(rbind(
  compute_fdr(truth = expected_label, test = df_var1$label, bool = FALSE),
  compute_fdr(truth = expected_label, test = df_var2$label, bool = FALSE),
  compute_fdr(truth = expected_label, test = df_var3$label, bool = FALSE)
  ))
colnames(df_rate) <- c("TP", "FN", "FP", "TN", "FDR", "POWER")
df_rate$model <- 1:3
write_tsv(df_rate, file = file.path(sdir, "beta_rate.tsv"))
rm(expected_label, df_rate)
##### 4.4 rank fdr and rank sensitivity #####
expected_label <- df_var1$expected_label
expected_label[expected_label == "ctrl"] <- "Neutral"
df_rankfdr <- data.frame(rbind(
  compute_rankfdr(true_label = expected_label, 
                  test_stats = df_var1$lfsr, 
                  test_value = df_var1$mean,
                  resolution = 100) %>% mutate(model = 1),
  compute_rankfdr(true_label = expected_label, 
                  test_stats = df_var2$lfsr, 
                  test_value = df_var2$mean,
                  resolution = 100) %>% mutate(model = 2),
  compute_rankfdr(true_label = expected_label, 
                  test_stats = df_var3$lfsr, 
                  test_value = df_var3$mean,
                  resolution = 100) %>% mutate(model = 3)
))
write_tsv(df_rankfdr, file = file.path(sdir, "beta_rankfdr.tsv"))
rm(expected_label)
# rankfdr plot
df_rankfdr_plot <- df_rankfdr %>% 
  select(model, rank, FDR, Sensitivity) %>%
  pivot_longer(cols = c(FDR, Sensitivity), names_to = "stats", values_to = "value")
p <- ggplot(df_rankfdr_plot, aes(x = rank, y = value)) +
  geom_step(aes(color = as.factor(model)), linewidth = 0.5, alpha = 0.5) +
  facet_wrap(~stats, scales = "free_y") +
  labs(color = NULL)
save_plot(file = file.path(sdir, "beta_rankfdr.png"), p, base_width = 10, base_height = 4)
rm(df_rankfdr_plot, df_rankfdr, p)
##### 4.5 fdr versus sensitivity #####
expected_label <- df_var1$expected_label
expected_label[expected_label == "ctrl"] <- "Neutral"
df_fdrsense <- data.frame(rbind(
  compute_fdrsense(true_label = expected_label, 
                   test_stats = df_var1$lfsr) %>% mutate(model = 1),
  compute_fdrsense(true_label = expected_label, 
                   test_stats = df_var2$lfsr) %>% mutate(model = 2),
  compute_fdrsense(true_label = expected_label, 
                   test_stats = df_var3$lfsr) %>% mutate(model = 3)
))
fdrsense_sig <- df_fdrsense %>% filter(fdr %in% c(0, 0.01, 0.05, 0.1))
write_tsv(df_fdrsense, file = file.path(sdir, "beta_fdrsense.tsv"))
p <- ggplot(df_fdrsense, aes(x = FDR, y = Power)) +
  geom_path(aes(color = as.factor(model)), alpha = 0.5) +
  labs(x = "false discovery rate", y = "power", color = NULL, shape = NULL) +
  geom_point(aes(x = FDR, y = Power, color = as.factor(model), shape = as.factor(fdr)), 
             data = fdrsense_sig, size = 2, alpha = 0.8) 
save_plot(file = file.path(sdir, "beta_fdrsense.png"), p, base_width = 5.5, base_height = 4)
rm(df_fdrsense, fdrsense_sig, p, expected_label)

########## 5 phi value ##########
##### expected phi #####
expected_phi <- (df_var1 %>% 
  group_by(pos, ctrl) %>% 
  summarise(expected_phi = mean(expected_effect)) %>%
  filter(ctrl == FALSE))[["expected_phi"]]
empr_phi1 <- (df_var1 %>% 
                group_by(pos, ctrl) %>% 
                summarise(empr_phi1 = mean(mean)) %>%
                filter(ctrl == FALSE))[["empr_phi1"]]
empr_phi2 <- (df_var2 %>% 
                group_by(pos, ctrl) %>% 
                summarise(empr_phi2 = mean(mean)) %>%
                filter(ctrl == FALSE))[["empr_phi2"]]
empr_phi3 <- (df_var3 %>% 
                group_by(pos, ctrl) %>% 
                summarise(empr_phi3 = mean(mean)) %>%
                filter(ctrl == FALSE))[["empr_phi3"]]

##### phi correlation #####
df_corr <- 
  data.frame(
    model = 1:3,
    pearson = c(
      cor(expected_phi, df_pos1$phi_mean, method = "pearson"),
      cor(expected_phi, df_pos2$phi_mean, method = "pearson"),
      cor(expected_phi, df_pos3$phi_mean, method = "pearson")
    ),
    spearman = c(
      cor(expected_phi, df_pos1$phi_mean, method = "spearman"),
      cor(expected_phi, df_pos2$phi_mean, method = "spearman"),
      cor(expected_phi, df_pos3$phi_mean, method = "spearman")
    )
  )
write_tsv(df_corr, file = file.path(sdir, "phi_corr.tsv"))
rm(df_corr)
##### phi correlation plot #####
df_phi <- data.frame(
  phi_model1 = df_pos1$phi_mean,
  phi_model2 = df_pos2$phi_mean,
  phi_model3 = df_pos3$phi_mean,
  phi_expected = expected_phi,
  phi_empr1 = empr_phi1,
  phi_empr2 = empr_phi2,
  phi_empr3 = empr_phi3
)
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
# p
save_plot(filename = file.path(sdir, "phi_corr_models.png"), 
          plot = p, base_height = 5, base_width = 15)
p1 <- ggplot(df_phi, aes(phi_expected, phi_model1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p2 <- ggplot(df_phi, aes(phi_expected, phi_model2)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p3 <- ggplot(df_phi, aes(phi_expected, phi_model3)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed")
p <- plot_grid(p1, p2, p3, nrow = 1)
rm(p1, p2, p3)
# p
save_plot(filename = file.path(sdir, "phi_corr_expected.png"), 
          plot = p, base_height = 5, base_width = 15)
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
# p
save_plot(filename = file.path(sdir, "phi_corr_empr.png"), 
          plot = p, base_height = 5, base_width = 15)
rm(df_phi, p, expected_phi)

########## 6 position-level estimate ##########
##### phi and sigma2 #####
p1 <- ggplot(df_pos1, aes(phi_mean, sigma2_mean)) +
  geom_point()
p2 <- ggplot(df_pos2, aes(phi_mean, sigma2_mean)) +
  geom_point()
p3 <- ggplot(df_pos3, aes(phi_mean, sigma2_mean)) +
  geom_point()
p <- plot_grid(p1, p2, p3, nrow = 1)
rm(p1, p2, p3)
save_plot(filename = file.path(sdir, "pos_phi_sigma2_scatter.png"), 
          plot = p, base_height = 5, base_width = 15)
##### model 3: phi and rho #####
p1 <- ggplot(df_pos3, aes(phi_mean, rho_mean)) +
  geom_point()
p2 <- ggplot(df_pos3, aes(phi_mean, sigma2_mean)) +
  geom_point()
p3 <- ggplot(df_pos3, aes(rho_mean, sigma2_mean)) +
  geom_point()
p <- plot_grid(p1, p2, p3, nrow = 1)
rm(p1, p2, p3)
save_plot(filename = file.path(sdir, "pos_model3_scatter.png"), 
          plot = p, base_height = 5, base_width = 15)
##### model 3: rho #####
if (sim_mode == 4) {
  df_rho <- df_var3 %>% filter(!ctrl) %>%
    select(pos, act4, rho_mean) %>% unique() 
  
  p1 <- ggplot(df_rho, aes(x = act4, y = rho_mean)) +
    geom_boxplot() +
    geom_jitter(alpha = 0.5, size = 0.1)
  p2 <- ggplot(df_rho, aes(rho_mean, fill = act4)) +
    geom_histogram(bins = 50, color = "grey")
  p <- plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 2))
  save_plot(filename = file.path(sdir, "rho_act_model3.png"), 
            plot = p, base_height = 5, base_width = 10)
  rm(p1, p2, p)
}

########## 7 nu value estimate ##########
########## nu value ##########
df_nu <- df_var3 %>% 
  select(blosum_score, nu_mean, nu_sd) %>% unique() %>% 
  arrange(blosum_score) %>%
  mutate(n = as.numeric(table(df_var3$blosum_score)))
write_tsv(df_nu, file = file.path(sdir, "df_nu.tsv"))

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
  labs(x = "blosum", y = "nu") +
  geom_text(aes(y = max(nu_mean), label = n), size = 2, vjust = -1)
p <- plot_grid(p1, p2, nrow = 1)
save_plot(filename = file.path(sdir, "nu_barplot.png"), 
          plot = p, base_height = 5, base_width = 13)
rm(p1, p2, p)

######### finale ##########
saveRDS(rosace, file = file.path(sdir, "rosace_eval_full.rds"))