library(tidyverse)
library(cowplot)

dir <- file.path("workspace", "domainase_0614")
ddir <- file.path(dir, "results_domainase", "analysis")
sdir <- file.path(dir, "results_domainase", "summary")

########## domain list ##########
domain_list <- read_table(file.path(dir, "data/domain_list.txt"), col_names = "dom_ID")
domain_list <- domain_list %>% 
  separate(dom_ID, c("gene", "family", "fpos"), remove = FALSE, sep = "_")  


########## protein family ##########
fam_count <- data.frame(table(domain_list$family))
colnames(fam_count) <- c("family", "count")
fam_count$family <- as.character(fam_count$family)
fam_count$fam_filter <- fam_count$family
fam_count <- fam_count %>%
  mutate(fam_filter = ifelse(count > 5, fam_filter, "others")) 

########## nu value PCA (cluster domain) ##########

##### read df_nu3 dataframe #####
df_nu3 <- data.frame()
for (i in 1:nrow(domain_list)) {
  key <- domain_list$dom_ID[i]
  df_nu3_sub <- read_tsv(file.path(ddir, key, "table", "df_nu3.tsv"))
  df_nu3 <- rbind(df_nu3, 
                  df_nu3_sub %>% 
                    mutate(dom_ID = key, family = domain_list$family[i]))
}
write_tsv(df_nu3, file = file.path(dir, "results_domainase", "summary", "df_nu3.tsv"))
rm(df_nu3_sub, ddir, key, i)

##### process df_nu3 #####
df_nu3 <- df_nu3 %>%
  rowwise() %>%
  mutate(lfsr.neg = stats::pnorm(0, mean = nu_mean, sd = nu_sd, lower.tail = FALSE),
         lfsr.pos = stats::pnorm(0, mean = nu_mean, sd = nu_sd, lower.tail = TRUE),
         lfsr = min(.data$lfsr.neg, .data$lfsr.pos)) %>%
  ungroup()
df_nu3 <- df_nu3 %>% left_join(fam_count %>% select(family, fam_filter))

##### widen df_nu3 #####
df_nu3_wide <- df_nu3 %>% 
  filter(!stop) %>%
  select(blosum_score, nu_mean, dom_ID, fam_filter) %>%
  pivot_wider(id_cols = c("dom_ID", "fam_filter"), names_from = "blosum_score", values_from = "nu_mean") %>%
  select(-`5`)
df_nu3_wide[is.na(df_nu3_wide)] <- 0
df_nu3_wide <- df_nu3_wide %>% rowwise() %>%
  mutate(nu_mean = mean(`-6`, `-5`, `-4`, `-3`, `-2`, `-1`, `0`, `1`, `2`, `3`),
         nu_sd = sd(c(`-6`, `-5`, `-4`, `-3`, `-2`, `-1`, `0`, `1`, `2`, `3`))) %>%
  ungroup()

##### standard deviation #####
p <- ggplot(df_nu3_wide, aes(nu_mean, nu_sd)) +
  geom_point(size = 0.8) +
  theme_classic()
save_plot(filename = file.path(sdir, "nu_sd_scatter.png"), 
          plot = p, base_height = 4, base_width = 5)
p <- ggplot(df_nu3_wide, aes(nu_sd)) +
  geom_histogram(aes(fill = fam_filter), color = "grey",
                 bins = 30, position = "identity") +
  theme_classic() +
  facet_wrap(vars(fam_filter)) 
save_plot(filename = file.path(sdir, "nu_sd_scatter_family.png"), 
          plot = p, base_height = 6.5, base_width = 10)
rm(p)

##### PCA #####
library(ggfortify)
pca_res <- prcomp(df_nu3_wide[3:12], center = FALSE, scale. = FALSE)
pca_res_center <- prcomp(df_nu3_wide[3:12], center = TRUE, scale. = FALSE)
pca_res_scale <- prcomp(df_nu3_wide[3:12], center = TRUE, scale. = TRUE)

### loading 
p1 <- autoplot(pca_res, data = df_nu3_wide, size = 0.5,
         loadings = TRUE, loadings.colour = 'red',
         loadings.label = TRUE, loadings.label.size = 3) +
  theme_classic() 
p2 <- autoplot(pca_res_center, data = df_nu3_wide, size = 0.5,
              loadings = TRUE, loadings.colour = 'red',
              loadings.label = TRUE, loadings.label.size = 3) +
  theme_classic() 
p3 <- autoplot(pca_res_scale, data = df_nu3_wide, size = 0.5,
               loadings = TRUE, loadings.colour = 'red',
               loadings.label = TRUE, loadings.label.size = 3) +
  theme_classic() 
p <- plot_grid(p1, p2, p3, nrow = 1)
save_plot(filename = file.path(sdir, "nu_mean_pca.png"), 
          plot = p, base_height = 4, base_width = 14)

### cluster by family
p <- autoplot(pca_res, data = df_nu3_wide, colour = 'fam_filter', size = 0.5) +
  theme_classic() +
  facet_wrap(vars(fam_filter))
save_plot(filename = file.path(sdir, "nu_mean_pca_family.png"), 
          plot = p, base_height = 6.5, base_width = 10)
p <- autoplot(pca_res_center, data = df_nu3_wide, colour = 'fam_filter', size = 0.5) +
  theme_classic() +
  facet_wrap(vars(fam_filter))
save_plot(filename = file.path(sdir, "nu_mean_pca_centered_family.png"), 
          plot = p, base_height = 6.5, base_width = 10)

rm(pca_res, pca_res_center, pca_res_scale)
rm(p1, p2, p3, p)

########## phi vs rho pattern for significant nu ##########
df_nu3  <- df_nu3 %>% rowwise() %>%
  mutate(test.sig = (lfsr <= 0.05),
         effect.sig = (abs(nu_mean) >= 0.5)) %>%
  ungroup()
domain_list$count_sig_test <- tapply(df_nu3$test.sig, df_nu3$dom_ID, sum)
domain_list$count_sig_effect <- tapply(df_nu3$effect.sig, df_nu3$dom_ID, sum)

table(domain_list$count_sig_test)
table(domain_list$count_sig_effect)




