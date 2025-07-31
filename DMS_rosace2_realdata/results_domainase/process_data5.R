library(tidyverse)
library(cowplot)

dir <- file.path("workspace", "domainase_0614")
ddir <- file.path(dir, "results_domainase", "analysis")
sdir <- file.path(dir, "results_domainase", "summary")

########## domain list ##########
domain_list <- read_table(file.path(dir, "data/domain_list.txt"), col_names = "dom_ID")
domain_list <- domain_list %>% 
  separate(dom_ID, c("gene", "family", "fpos"), remove = FALSE, sep = "_")  


##### load variant data frame and compute summary statistics #####
df_decomp <- data.frame()
for (i in 1:nrow(domain_list)) {
  key <- domain_list$dom_ID[i]
  load(file.path(ddir, key, "df_var_full.rda"))
  df_decomp_key <- data.frame(
    dom_ID = key,
    model = 1:3,
    var_phi = c(
      var(df_var1$phi_mean) / var(df_var1$mean),
      var(df_var2$phi_mean) / var(df_var2$mean),
      var(df_var3$phi_mean) / var(df_var3$mean)),
    var_nu = c(
      0,
      var(df_var2$nu_mean) / var(df_var2$mean),
      var(df_var3$nu_mean * df_var3$rho_mean) / var(df_var3$mean)
    )
  )
  df_decomp <- rbind(df_decomp, df_decomp_key)
  rm(df_var1, df_var2, df_var3)
}
df_decomp <- df_decomp %>% mutate(var_all = var_phi + var_nu)
write_tsv(df_decomp, file = file.path(dir, "results_domainase", "summary", "df_decomp.tsv"))
rm(i, key, df_decomp_key)

##### summary statistics #####
df_decomp <- read_tsv(file.path(dir, "results_domainase", "summary", "df_decomp.tsv"))

### choice 1: figure 2 ###
example <- c("O60341_PF04433_186", "P12081_PF00458_8", "P19793_PF00105_132", "P49792_PF00641_1478")
df_decomp %>% filter(dom_ID %in% example)

for (key in example) {
  df_decomp_sub <- df_decomp %>% filter(dom_ID == key) %>%
    select(-var_all) %>%
    pivot_longer(cols = 3:4)
  df_decomp_sub$model <- rep(c("position only", 
                               "w/ AA effect",
                               "w/ position-scaled AA effect"), each = 2)
  df_decomp_sub$name[4] <- "df_nu_no_scaling"
  df_decomp_sub$name <- factor(df_decomp_sub$name)
  levels(df_decomp_sub$name) <- c("unscaled AA effect", "position-scaled AA effect", "position effect")
  p <- ggplot(df_decomp_sub, aes(fill=name, y=value, x=as.factor(model))) + 
    geom_bar(position="stack", stat="identity", alpha = 0.8, color = "grey") + 
    coord_flip(ylim = c(0.5, 1)) +
    scale_fill_manual(values=c("#029E73", '#a369b0', '#0173B2')) +
    labs(x = "", y = "proportion of variance explained") +
    guides(fill = "none") +
    theme_classic()
  cowplot::save_plot(filename = file.path(dir, "results_domainase", "summary", "decomp", paste("decomp", key, ".png", sep = "")),
                     plot = p, base_height = 2, base_width = 7)
}

ggplot(df_decomp_sub, aes(fill=name, y=value, x=as.factor(model))) + 
  geom_bar(position="stack", stat="identity", alpha = 0.8, color = "grey") + 
  coord_flip(ylim = c(0.5, 1)) +
  scale_fill_manual(values=c("#029E73", '#a369b0', '#0173B2')) +
  labs(x = "", y = "proportion of variance explained") +
  theme_classic()
rm(p, df_decomp_sub, key)

### choice 2: global trend ###
df_decomp_wide <- df_decomp %>%
  pivot_wider(names_from = model, values_from = 3:5)

p <- ggplot(df_decomp_wide, aes(var_phi_1, var_phi_3)) +
  geom_point(color = '#0173B2') +
  geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 1, linetype = "dashed") +
  theme_classic() +
  labs(x = "variance explained by position effect (model 1)", 
       y = "variance explained by position effect (model 3)")
cowplot::save_plot(filename = file.path(dir, "results_domainase", "summary", "decomp", 
                                        paste("decomp", "position", ".png", sep = "")),
                   plot = p, base_height = 4, base_width = 5)


p <- ggplot(df_decomp_wide %>%
         mutate(color = ifelse(var_nu_3 > var_nu_2 + 0.01, "3", "0"),
                color = ifelse(var_nu_2 > var_nu_3 + 0.01, "2", color)), 
       aes(var_nu_2, var_nu_3)) +
  geom_point(aes(color = color)) +
  scale_color_manual(values=c('grey', "#029E73", '#a369b0')) +
  theme_classic() +
  labs(x = "variance explained by AA effect (model 2)", 
       y = "variance explained by AA effect (model 3)") +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  guides(color = "none")
cowplot::save_plot(filename = file.path(dir, "results_domainase", "summary", "decomp", 
                                        paste("decomp", "AA", ".png", sep = "")),
                   plot = p, base_height = 4, base_width = 5)

##### summary statistics -- analysis on protein family #####
df_decomp <- read_tsv(file.path(dir, "results_domainase", "summary", "df_decomp.tsv"))
df_decomp_family <- df_decomp %>% left_join(domain_list)

library(forcats)

df_decomp_family_plot <- df_decomp_family %>% 
  filter(model == 3) %>%
  group_by(family) %>% mutate(family_count = n()) %>% ungroup() %>%
  filter(family_count > 5)

df_decomp_family_plot_sum <- df_decomp_family_plot %>% group_by(family) %>% summarise(mean(var_phi), mean(var_nu))
df_decomp_family_plot$family <- fct_reorder(df_decomp_family_plot$family, 
                                            df_decomp_family_plot$var_phi, .fun = mean)
p1 <- ggplot(df_decomp_family_plot, aes(x = family, y = var_phi)) +
  geom_boxplot(outlier.shape = NA) +               # Box plot without showing outliers
  geom_jitter(width = 0.2, color = '#0173B2', alpha = 0.6) +  # Add individual points with some jitter
  labs(x = "Domain Family", y = "Variance explained by position effect (model 3)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p2 <- ggplot(df_decomp_family_plot, aes(x = family, y = var_nu)) +
  geom_boxplot(outlier.shape = NA) +               # Box plot without showing outliers
  geom_jitter(width = 0.2, color = '#a369b0', alpha = 0.6) +  # Add individual points with some jitter
  labs(x = "Domain Family", y = "Variance explained by AA effect (model 3)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
p <- plot_grid(p1, p2, nrow = 1)
cowplot::save_plot(filename = file.path(dir, "results_domainase", "summary", "decomp", 
                                        paste("decomp", "AA", "domain", ".png", sep = "")),
                   plot = p, base_height = 4, base_width = 14)



# p1 <- ggplot(df_decomp %>% filter(model == 1)) +
#   geom_histogram(aes(var_phi), fill = '#0173B2', alpha = 0.5) + 
#   theme_classic() +
#   labs(x = "variance explained by position effect (model 1)")
# p2 <- ggplot(df_decomp %>% 
#          filter(model %in% c(2, 3))) +
#   geom_histogram(aes(var_nu, fill = factor(model)), position = 'identity',
#                  alpha = 0.5) +
#   scale_fill_manual(values=c("#029E73", '#a369b0')) +
#   theme_classic() +
#   labs(x = "variance explained by AA effect (model 2 or 3)",
#        fill = "") +
#   guides(fill = "none")
# plot_grid(p1, p2, nrow = 2)



# p2 <- ggplot(df_decomp %>% filter(model == 2, m)) +
#   geom_histogram(aes(var_nu), fill = "#029E73", color = "white") +
#   theme_classic() +
#   labs(x = "variance explained by AA effect (model 2)")
# p3 <- ggplot(df_decomp %>% filter(model == 3)) +
#   geom_histogram(aes(var_nu), fill = '#a369b0', color = "white") + 
#   theme_classic() +
#   labs(x = "variance explained by position-scaled AA effect (model 3)")
# 






