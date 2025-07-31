# library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(cowplot)

base <- "workspace/downsample_0817"
dir <- "workspace/downsample_0817/data"

load(file.path(dir, "oct1_sm73.rda"))
rosace_full <- rosace

########### random dataset ###########
size_list <- c(80, 60, 40, 20)
for (i in 1:5) {
  set.seed(i)
  for (size in size_list) {
    selected_vars <- rosace_full@var.data %>%
      group_by(pos) %>%
      sample_frac(size=size/100) %>%
      select(variants)
    selected_rows <- sort(match(selected_vars$variants, rosace_full@assay.sets$`1SM73`@var.names))
    rosace <- rosace_full
    rosace@assay.sets$`1SM73`@raw.counts <- rosace@assay.sets$`1SM73`@raw.counts[selected_rows, ]
    rosace@assay.sets$`1SM73`@combined.counts <- rosace@assay.sets$`1SM73`@combined.counts[selected_rows, ]
    rosace@assay.sets$`1SM73`@var.names <- rosace@assay.sets$`1SM73`@var.names[selected_rows]
    save(rosace, file = file.path(dir, paste("oct1_sm73_random", size, "_", i, ".rda", sep = "")))
  }
}
rm(i, size, selected_vars, selected_rows, rosace)

########### adjusted dataset ###########
# load(file.path("workspace/downsample_0817", "oct1_sm73", "df_var_full.rda"))
# df_var <- df_var3 
# rm(df_var1, df_var2, df_var3)
# df_var_start <- df_var %>% select(1:7, blosum_score)

##### choice round 1 #####
# for (i in 1:5) {
#   set.seed(i)
#   selected_vars <- rbind(
#     df_var_start %>% 
#       filter(blosum_score %in% c(5, -7)) %>%
#       sample_frac(size=0.4),
#     df_var_start %>% 
#       filter(blosum_score != 5, blosum_score != -7) %>%
#       group_by(pos, blosum_score) %>%
#       sample_n(size = 1)
#   )
#   # nrow(selected_vars)/nrow(df_var_start) # 0.36
#   selected_rows <- sort(match(selected_vars$variants, rosace_full@assay.sets$`1SM73`@var.names))
#   rosace <- rosace_full
#   rosace@assay.sets$`1SM73`@raw.counts <- rosace@assay.sets$`1SM73`@raw.counts[selected_rows, ]
#   rosace@assay.sets$`1SM73`@combined.counts <- rosace@assay.sets$`1SM73`@combined.counts[selected_rows, ]
#   rosace@assay.sets$`1SM73`@var.names <- rosace@assay.sets$`1SM73`@var.names[selected_rows]
#   save(rosace, file = file.path(dir, paste("oct1_sm73_choice1", "_", i, ".rda", sep = "")))
# }
# rm(selected_rows, selected_vars, rosace, i)

##### choice round 2 #####
# strategy: 
# high amino acid sensitive position: sample more? how much
# low amino acid sensitive position: sample less? how less
# weights
# ddir <- file.path(base, "results_oct", "analysis", "oct1_sm73_choice1_1")


########## analyze full data ########## 

###### nu plot
df_nu <- read_tsv(file.path("workspace/downsample_0817", 
                            "oct1_sm73", "table", "df_nu3.tsv"))
df_nu <- df_nu %>% filter(blosum_score != -7, blosum_score != 5)

p <- ggplot(df_nu, aes(x = as.factor(blosum_score))) +
  geom_bar(aes(y = nu_mean), stat = "identity", 
           fill = "#029E73", alpha = 0.7) +
  geom_point(aes(y = nu_mean), colour="darkgreen", alpha = 0.9) +
  geom_errorbar(aes(ymin = nu_mean - nu_sd, ymax = nu_mean + nu_sd), 
                width=0.3, colour="darkgreen", alpha=0.9, size=1) +
  labs(x = "blosum", y = "nu") +
  # geom_text(aes(y = max(nu_mean), label = n), size = 2, vjust = -1) +
 # coord_flip()+
  theme_classic()
cowplot::save_plot(filename = file.path(base, "plot", "main_plot4_nu.png"), 
                   plot = p, base_height = 3, base_width = 6)

######
load(file.path("workspace/downsample_0817", "oct1_sm73", "df_pos_full.rda"))
# df_pos <- df_pos3
# rm(df_pos1, df_pos2, df_pos3)
df_sigma2 <- data.frame(
  model1_sigma2 = df_pos1$sigma2_mean,
  model2_sigma2 = df_pos2$sigma2_mean,
  model3_sigma2 = df_pos3$sigma2_mean
)
df_sigma2 <- df_sigma2 %>%
  mutate(label_sigma2_diff = (model1_sigma2 - model3_sigma2) > 0.1)
p <- ggplot(df_sigma2, aes(model1_sigma2, model3_sigma2)) +
  geom_point(aes(color = label_sigma2_diff)) +
  theme_classic() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "sigma2 (model1)", y = "sigma2 (model3)") +
  guides(color = "none") +
  scale_color_manual(values=c('#999999','#a369b0'))
cowplot::save_plot(filename = file.path(base, "plot", "main_plot1_sigma2.png"), 
                   plot = p, base_height = 4, base_width = 5)

df_sigma2 <- df_sigma2 %>%
  mutate(diff_sigma2_model1_vs_3 = model1_sigma2 - model3_sigma2,
         model3_rho = df_pos3$rho_mean) 
p <- ggplot(df_sigma2, aes(diff_sigma2_model1_vs_3, model3_rho)) +
  geom_point(aes(color = label_sigma2_diff)) +
  theme_classic() +
  labs(x = "sigma2 (model1 - model3)", y = "rho (model3)") +
  guides(color = "none") +
  scale_color_manual(values=c('#999999','#a369b0'))
cowplot::save_plot(filename = file.path(base, "plot", "main_plot1_sigma2_rho.png"), 
                   plot = p, base_height = 4, base_width = 5)


# df_pos3 <- df_pos3 %>%
#   mutate(label_not_blosum_AA = ifelse(sigma2_mean > 0.8, TRUE, FALSE))
p <- ggplot(df_pos3, aes(rho_mean, sigma2_mean)) +
  geom_point(aes(color = sigma2_mean)) +
  theme_classic() +
  scale_color_gradient(low = "grey", high = "darkorange") +
  labs(x = "rho (model 3)", y = "sigma2 (model 3)") +
  guides(color = "none")
cowplot::save_plot(filename = file.path(base, "plot", "main_plot1_sigma2_rho2.png"), 
                   plot = p, base_height = 4, base_width = 5)

##########  full data tier ##########  
load(file.path("workspace/downsample_0817", "oct1_sm73", "df_var_full.rda"))
df_var <- df_var3
rm(df_var1, df_var2, df_var3)
df_var <- df_var %>% mutate(pred = phi_mean + rho_mean * nu_mean)

# full data: position tier dropping 
df_var_missense <- df_var %>% filter(type == "missense")
breaks <- quantile(df_var_missense$rho_mean, probs = seq(0, 1, by = 0.05), na.rm = TRUE)
df_var_missense$rho_group <- cut(df_var_missense$rho_mean, breaks = breaks, include.lowest = TRUE, labels = FALSE)

breaks <- quantile(df_var_missense$sigma2_mean, probs = seq(0, 1, by = 0.05), na.rm = TRUE)
df_var_missense$sigma2_group <- cut(df_var_missense$sigma2_mean, breaks = breaks, include.lowest = TRUE, labels = FALSE)

breaks <- quantile(df_var_missense$phi_mean, probs = seq(0, 1, by = 0.05), na.rm = TRUE)
df_var_missense$phi_group <- cut(df_var_missense$phi_mean, breaks = breaks, include.lowest = TRUE, labels = FALSE)

df_rsquared_rho_tier <- data.frame()
for (i in 1:10) {
  fit1 <- summary(lm(mean ~ pred, data = df_var_missense %>% filter(rho_group == i)))
  fit2 <- summary(lm(mean ~ pred, data = df_var_missense %>% filter(sigma2_group == i)))
  temp <- data.frame(
    type = c("rho", "sigma2"),
    tier = c(i/10, i/10),
    r.squared = c(fit1$r.squared, fit2$r.squared),
    adj.r.squared = c(fit1$adj.r.squared, fit2$adj.r.squared)
  )
  df_rsquared_rho_tier <- rbind(df_rsquared_rho_tier, temp)
}
rm(fit1, fit2, temp, i)

p <- ggplot(df_rsquared_rho_tier, aes(as.factor(tier), r.squared)) +
  geom_col(aes(fill = type), position = "dodge", alpha = 0.5) +
  theme_classic() +
  facet_wrap(vars(type), nrow = 2) +
  guides(fill = "none") + 
  ylim(0, 1) +
  labs(x = "percentile (low to high)") +
  scale_fill_manual(values=c('#a369b0', 'darkorange')) +
  geom_hline(yintercept = 0.6510702, color = "black", linetype = "dashed")
cowplot::save_plot(filename = file.path(base, "plot", "main_plot2_rsquared.png"), 
                   plot = p, base_height = 4, base_width = 5)
rm(p)

########### compare to alphamissense ###########
df_alpham <- read_tsv(file.path(base, "oct1_alpham.tsv"))
df_alpham <-df_alpham %>% 
  select(pos = position, wt = a.a.1, mut = a.a.2, 
         alpham = `pathogenicity score`, alpham_class = `pathogenicity class`)
df_var_missense <- df_var_missense %>% left_join(df_alpham)

df_rsquared_alpham_tier <- data.frame()
for (i in 1:10) {
  fit1 <- summary(lm(mean ~ alpham, data = df_var_missense %>% filter(rho_group == i)))
  fit2 <- summary(lm(mean ~ alpham, data = df_var_missense %>% filter(sigma2_group == i)))
  fit3 <- summary(lm(mean ~ alpham, data = df_var_missense %>% filter(phi_group == i)))
  temp <- data.frame(
    type = c("rho", "sigma2", "phi"),
    tier = c(i/10, i/10, i/10),
    r.squared = c(fit1$r.squared, fit2$r.squared, fit3$r.squared),
    adj.r.squared = c(fit1$adj.r.squared, fit2$adj.r.squared, fit3$adj.r.squared)
  )
  df_rsquared_alpham_tier <- rbind(df_rsquared_alpham_tier, temp)
}
rm(fit1, fit2, temp, i)

p <- ggplot(df_rsquared_alpham_tier %>% filter(type != "phi"), aes(as.factor(tier), r.squared)) +
  geom_col(aes(fill = type), position = "dodge", alpha = 0.5) +
  theme_classic() +
  facet_wrap(vars(type), nrow = 3) +
  guides(fill = "none") +
  labs(x = "percentile (low to high)") +
  ylim(0, 1) +
  scale_fill_manual(values=c('#a369b0', 'darkorange')) +
  geom_hline(yintercept = 0.6510702, color = "black", linetype = "dashed")
cowplot::save_plot(filename = file.path(base, "plot", "supp_plot2_rsquared.png"), 
                   plot = p, base_height = 4, base_width = 5)

########## correlation analysis -- review added ########## 
# Compute correlation score between mean and alpha_m for each position
df_correlation_rho <- df_var_missense %>%
  group_by(rho_group) %>%
  summarize(correlation_alpham = cor(mean, alpham, use = "complete.obs")^2,
            correlation_rosace = cor(mean, pred, use = "complete.obs")^2,
            phi_mean = mean(phi_mean, na.rm = TRUE),
            rho_mean = mean(rho_mean, na.rm = TRUE),
            sigma2_mean = mean(sigma2_mean, na.rm = TRUE),
            alpham_perc_ambiguous = mean(alpham_class == "ambiguous"),
            alpham_perc_benign = mean(alpham_class == "likely_benign"),
            alpham_perc_patho = mean(alpham_class == "likely_pathogenic"),
            ) %>%
  ungroup()

df_correlation_sigma2 <- df_var_missense %>%
  group_by(sigma2_group) %>%
  summarize(correlation_alpham = cor(mean, alpham, use = "complete.obs")^2,
            correlation_rosace = cor(mean, pred, use = "complete.obs")^2,
            phi_mean = mean(phi_mean, na.rm = TRUE),
            rho_mean = mean(rho_mean, na.rm = TRUE),
            sigma2_mean = mean(sigma2_mean, na.rm = TRUE),
            alpham_perc_ambiguous = mean(alpham_class == "ambiguous"),
            alpham_perc_benign = mean(alpham_class == "likely_benign"),
            alpham_perc_patho = mean(alpham_class == "likely_pathogenic"),
  ) %>%
  ungroup()

write_csv(df_correlation_rho, "correlation_rho.csv")
write_csv(df_correlation_sigma2, "correlation_sigma2.csv")

summary(lm(correlation_alpham ~ sigma2_mean, df_correlation_sigma2)) # <2e-16, 0.1888
summary(lm(correlation_rosace ~ sigma2_mean, df_correlation_sigma2)) # 0.0708, 0.005935

summary(lm(correlation_alpham ~ rho_mean, df_correlation_rho))  # <2e-16, 0.4472
summary(lm(correlation_rosace ~ rho_mean, df_correlation_rho)) # <2e-16, 0.2134

# Plot the correlation between phi_mean and the computed correlation
p1 <- ggplot(df_correlation_sigma2, aes(x = sigma2_mean, y = correlation_alpham)) +
  geom_point(alpha = 1, color = 'darkorange') +
  theme_classic() +
  labs(x = "Sigma2 Mean", y = "R2 (AlphaMissense)") 
  # geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE)

p2 <- ggplot(df_correlation_rho, aes(x = rho_mean, y = correlation_alpham)) +
  geom_point(alpha = 1, color = '#a369b0') +
  theme_classic() +
  labs(x = "Rho Mean", y = "R2 (AlphaMissense)") 
  # geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE)

p3 <- ggplot(df_correlation_sigma2, aes(x = sigma2_mean, y = correlation_rosace)) +
  geom_point(alpha = 1, color = 'darkorange') +
  theme_classic() +
  labs(x = "Sigma2 Mean", y = "R2 (Posterior)") 
  # geom_smooth(method = "lm", color = "red", linetype = "dashed", se = FALSE)

p4 <- ggplot(df_correlation_rho, aes(x = rho_mean, y = correlation_rosace)) +
  geom_point(alpha = 1, color = '#a369b0') +
  theme_classic() +
  labs(x = "Rho Mean", y = "R2 (Posterior)") 

# cowplot
plot_grid(p4, p2, p3, p1)

# 
# 
ggplot(df_correlation_sigma2, aes(x = correlation_alpham, y = alpham_perc_ambiguous)) +
  geom_point(alpha = 0.7, color = "darkblue") +
  theme_classic() +
  labs(x = "Sigma2 Mean", y = "R2 (DMS Score vs AlphaMissense Prediction)",
       title = "Correlation between Sigma2 and Alphamissense R2") 
summary(lm(alpham_perc_ambiguous~correlation_alpham, df_correlation_sigma2))


# cowplot::save_plot(filename = file.path(base, "plot", "correlation_phi_alpha_m.png"), 
#                    plot = p, base_height = 4, base_width = 6)



########## process random 20 data ########## 

for (i in 1:5) {
  df_var_full <- read_tsv(
    file.path(base, "results_oct", "analysis", 
              paste("oct1_sm73_random20", i, sep = "_"), 
              "table", "df_var_full.tsv")
  )
  df_var_full_missing <- df_var_full %>% filter(is.na(mean_sample))
  breaks <- quantile(df_var_full_missing$sigma2_mean, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  df_var_full_missing$sigma2_group <- cut(df_var_full_missing$sigma2_mean, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  # nrow(df_var_full) * 0.2 # 2075
  set.seed(i)
  selected_vars <- rbind(
    df_var_full_missing %>% filter(sigma2_group == 1) %>% group_by(pos) %>% sample_frac(2/20),
    df_var_full_missing %>% filter(sigma2_group == 2) %>% group_by(pos) %>% sample_frac(4/20),
    df_var_full_missing %>% filter(sigma2_group == 3) %>% group_by(pos) %>% sample_frac(6/20),
    df_var_full_missing %>% filter(sigma2_group == 4) %>% group_by(pos) %>% sample_frac(8/20)
  ) %>% ungroup() %>% select(variants)
  selected_vars <- rbind(
    selected_vars,
    df_var_full %>% filter(!is.na(mean_sample)) %>% select(variants)
  )
  selected_rows <- sort(match(selected_vars$variants, rosace_full@assay.sets$`1SM73`@var.names))
  rosace <- rosace_full
  rosace@assay.sets$`1SM73`@raw.counts <- rosace@assay.sets$`1SM73`@raw.counts[selected_rows, ]
  rosace@assay.sets$`1SM73`@combined.counts <- rosace@assay.sets$`1SM73`@combined.counts[selected_rows, ]
  rosace@assay.sets$`1SM73`@var.names <- rosace@assay.sets$`1SM73`@var.names[selected_rows]
  save(rosace, file = file.path(dir, paste("oct1_sm73_choice", 40, "_", i, ".rda", sep = "")))
}

########## process choice 40 data ########## 

for (i in 1:5) {
  df_var_full <- read_tsv(
    file.path(base, "results_oct", "analysis", 
              paste("oct1_sm73_choice40", i, sep = "_"), 
              "table", "df_var_full.tsv")
  )
  df_var_full_missing <- df_var_full %>% filter(is.na(mean_sample))
  breaks <- quantile(df_var_full_missing$sigma2_mean, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  df_var_full_missing$sigma2_group <- cut(df_var_full_missing$sigma2_mean, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  # nrow(df_var_full) * 0.6 # 2075
  set.seed(i)
  selected_vars <- rbind(
    df_var_full_missing %>% filter(sigma2_group == 1) %>% group_by(pos) %>% sample_n(2, replace = TRUE),
    df_var_full_missing %>% filter(sigma2_group == 2) %>% group_by(pos) %>% sample_n(3, replace = TRUE),
    df_var_full_missing %>% filter(sigma2_group == 3) %>% group_by(pos) %>% sample_n(6, replace = TRUE),
    df_var_full_missing %>% filter(sigma2_group == 4) %>% group_by(pos) %>% sample_n(7, replace = TRUE)
  ) %>% ungroup() %>% select(variants) %>% unique()
  selected_vars <- rbind(
    selected_vars,
    df_var_full %>% filter(!is.na(mean_sample)) %>% select(variants)
  )
  print(nrow(selected_vars))
  selected_rows <- sort(match(selected_vars$variants, rosace_full@assay.sets$`1SM73`@var.names))
  rosace <- rosace_full
  rosace@assay.sets$`1SM73`@raw.counts <- rosace@assay.sets$`1SM73`@raw.counts[selected_rows, ]
  rosace@assay.sets$`1SM73`@combined.counts <- rosace@assay.sets$`1SM73`@combined.counts[selected_rows, ]
  rosace@assay.sets$`1SM73`@var.names <- rosace@assay.sets$`1SM73`@var.names[selected_rows]
  save(rosace, file = file.path(dir, paste("oct1_sm73_choice", 60, "_", i, ".rda", sep = "")))
}

########## process choice 60 data ########## 

for (i in 1:5) {
  df_var_full <- read_tsv(
    file.path(base, "results_oct", "analysis", 
              paste("oct1_sm73_choice60", i, sep = "_"), 
              "table", "df_var_full.tsv")
  )
  df_var_full_missing <- df_var_full %>% filter(is.na(mean_sample))
  breaks <- quantile(df_var_full_missing$sigma2_mean, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
  df_var_full_missing$sigma2_group <- cut(df_var_full_missing$sigma2_mean, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  
  # nrow(df_var_full) * 0.8 # 8302
  set.seed(i)
  selected_vars <- rbind(
    df_var_full_missing %>% filter(sigma2_group == 1) %>% group_by(pos) %>% sample_n(2, replace = TRUE),
    df_var_full_missing %>% filter(sigma2_group == 2) %>% group_by(pos) %>% sample_n(5, replace = TRUE),
    df_var_full_missing %>% filter(sigma2_group == 3) %>% group_by(pos) %>% sample_n(7, replace = TRUE),
    df_var_full_missing %>% filter(sigma2_group == 4) %>% group_by(pos) %>% sample_n(10, replace = TRUE)
  ) %>% ungroup() %>% select(variants) %>% unique()
  selected_vars <- rbind(
    selected_vars,
    df_var_full %>% filter(!is.na(mean_sample)) %>% select(variants)
  )
  print(nrow(selected_vars))
  selected_rows <- sort(match(selected_vars$variants, rosace_full@assay.sets$`1SM73`@var.names))
  rosace <- rosace_full
  rosace@assay.sets$`1SM73`@raw.counts <- rosace@assay.sets$`1SM73`@raw.counts[selected_rows, ]
  rosace@assay.sets$`1SM73`@combined.counts <- rosace@assay.sets$`1SM73`@combined.counts[selected_rows, ]
  rosace@assay.sets$`1SM73`@var.names <- rosace@assay.sets$`1SM73`@var.names[selected_rows]
  save(rosace, file = file.path(dir, paste("oct1_sm73_choice", 80, "_", i, ".rda", sep = "")))
}


########### analyze random dataset ###########

##### df_rsquared #####
df_rsquared <- read_tsv(file.path(base, "results_oct", "summary", "df_rsquared.tsv"))
df_rsquared <- df_rsquared %>% 
  separate(run_ID, into = c("gene", "drug", "sample", "rep")) %>% 
  select(-gene, -drug, -rep) %>%
  filter(sample != "choice1")

##### add full data row #####
load(file.path("workspace/downsample_0817", "oct1_sm73", "df_var_full.rda"))
df_var <- df_var3
rm(df_var1, df_var2, df_var3)
df_var <- df_var %>% mutate(pred = phi_mean + rho_mean * nu_mean)

# fit model r squared
fit <- summary(lm(mean ~ pred, data = df_var)) # 0.6511
temp <- data.frame(fit$r.squared, fit$adj.r.squared, 0, 0, "random100")
colnames(temp) <- colnames(df_rsquared)
df_rsquared <- rbind(df_rsquared, temp)
rm(fit, temp)

### plot 1
df_rsquared_random_plot <- df_rsquared %>% 
  filter(substr(sample, 1, 6) == "random") %>%
  mutate(perc = as.numeric(substr(sample, 7, nchar(sample)))/100) %>%
  mutate(key = "random")
df_rsquared_choice_plot <- df_rsquared %>% 
  filter(substr(sample, 1, 6) == "choice") %>%
  mutate(perc = as.numeric(substr(sample, 7, nchar(sample)))/100) %>%
  mutate(key = "sigma2-guided")
df_rsquared_plot <- rbind(df_rsquared_random_plot, df_rsquared_choice_plot)
rm(df_rsquared_random_plot, df_rsquared_choice_plot)


p <- ggplot(df_rsquared_plot %>% filter(sample != "random100"), 
       aes(as.factor(perc), r.squared_missing)) +
  geom_jitter(aes(color = key)) +
  geom_boxplot(aes(color = key), position = "identity") +
  geom_hline(yintercept = 0.6510702, linetype = "dashed") +
  ylim(0, 1) +
  theme_classic() +
  scale_color_manual(values=c('grey', 'darkorange')) +
  labs(color = "", x = "percentage of variants learned",
       y = "R-sqaured of variants missing")
cowplot::save_plot(filename = file.path(base, "plot", "main_plot3_rsquared.png"), 
                   plot = p, base_height = 4, base_width = 6)








