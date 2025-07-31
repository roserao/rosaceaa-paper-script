library(tidyverse)
library(cowplot)

dir <- "workspace/met_0812/results_met/analysis"
sdir_plot <- "workspace/met_0812/plot_0815"
cond_list <- list.dirs(dir, full.names = FALSE, recursive = FALSE)

############ position full ##########
df_pos_full <- data.frame()
for (cond in cond_list) {
  load(file.path(dir, cond, "df_pos_full.rda"))
  df_pos_full <- rbind(df_pos_full, df_pos3 %>% mutate(cond = cond))
  rm(df_pos1, df_pos2, df_pos3)
}

df_pos_full <- df_pos_full %>% separate(cond, into = c("gene", 'inh', "genetic"), remove = FALSE)
df_pos_full$pos <- df_pos_full$pos + 1058

############ position full ##########
criz_pocket <- c(1108, 1159, 1163, 1260, 1158, 1211, 1092, 1230)
criz_pocket <- sort(criz_pocket)
atp_pocket <- c(1089, 1159, 1260, 1222, 1223, 1227)

############ by condition ##########
ggplot(df_pos_full, aes(phi_mean, rho_mean)) +
  geom_point(size = 0.5) +
  facet_wrap(vars(cond), ncol = 6) +
  theme_classic()

ggplot(df_pos_full %>% 
         filter(inh == "Crizo") %>%
         mutate(crizo_pocket = pos %in% criz_pocket), 
       aes(phi_mean, rho_mean)) +
  geom_point(aes(color = crizo_pocket), size = 0.5) +
  facet_wrap(vars(cond), ncol = 2) +
  theme_classic()

############ main plot 1: clustering ##########
df_pos_summarise <- df_pos_full %>% 
  filter(genetic == "WT") %>%
  group_by(pos) %>%
  summarise(phi_mean_cond = mean(phi_mean), phi_var_cond = var(phi_mean),
            rho_mean_cond = mean(rho_mean), rho_var_cond = var(rho_mean))
df_pos_summarise <- df_pos_summarise %>%
  mutate(cluster = ifelse(rho_var_cond < 0.01 & phi_var_cond < 1, 0, -1),
         cluster = ifelse(rho_var_cond >= 0.01 & phi_var_cond >= 1, 3, cluster),
         cluster = ifelse(rho_var_cond >= 0.01 & phi_var_cond < 1, 2, cluster),
         cluster = ifelse(rho_var_cond < 0.01 & phi_var_cond >= 1, 1, cluster))
table(df_pos_summarise$cluster)

p <- ggplot(df_pos_summarise, aes(phi_var_cond, rho_var_cond)) +
  geom_point(aes(color = as.factor(cluster))) +
  geom_hline(yintercept = 0.01, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  scale_color_manual(values=c('#999999','#0173B2', '#a369b0', "#ffc48a")) +
  labs(x = "variance of position score (phi) across inhibitors", y = "variance of AA scaling (rho) across inhibitors") +
  guides(fill = "none", color = "none", linetype = "none", shape = "none") +
  theme_classic()
cowplot::save_plot(filename = file.path(sdir_plot, "main_plot1_cluster.png"), 
                   plot = p, base_height = 4, base_width = 5)

############ main plot 2: position ##########
pos_grey <- c(1327)
pos_purple <- c(1290, 1188, 1092, 1326, 1219, 1113, 1063)
pos_orange <- c(1314, 1294, 1284)
pos_blue <- c(1163, 1158, 1230, 1175, 1166, 1291, 1249)

sort((df_pos_summarise %>% filter(cluster == 2))$pos)

df_pos_full[df_pos_full$inh == "Camp", "inh"] <- "Cap"
df_inh <- data.frame(
  inh = c("A458", "Cabo", "Cap", "Crizo", "DMSO", "Gle", "Glu", "Mere", "NVP", "Savo", "Tepo", "Tiv"),
  inh_type = c("I1/2", "II", "Ia", "Ia", "DMSO", "II", "Ia", "II", "Ib", "Ib", "Ia", "III")
)
df_pos_full <- df_pos_full %>% left_join(df_inh)

color_grey <- '#999999'
color_orange <- "#ffc48a"
color_purple <- '#a369b0'
color_blue <- '#0173B2'

df_plot_pos <- data.frame(
  pos_idx = c(pos_grey, pos_purple, pos_orange, pos_blue),
  color_idx = c(rep(color_grey, length(pos_grey)), 
                rep(color_purple, length(pos_purple)), 
                rep(color_orange, length(pos_orange)), 
                rep(color_blue, length(pos_blue)))
)

for (i in 1:nrow(df_plot_pos)) {
  color_idx <- df_plot_pos$color_idx[i]
  pos_idx <- df_plot_pos$pos_idx[i]
  p <- ggplot(df_pos_full %>% filter(pos == pos_idx, genetic == "WT") %>%
                mutate(ctrl = inh == "DMSO"), aes(phi_mean, rho_mean)) +
    # geom_point(aes(shape = ctrl), colour = color_idx, size = 2) +
    geom_point(aes(shape = ctrl, colour = inh_type), size = 2) +
    ggrepel::geom_text_repel(aes(label = inh)) +
    guides(shape = "none", color = "none") +
    xlim(min(df_pos_full$phi_mean), max(df_pos_full$phi_mean)) +
    ylim(min(df_pos_full$rho_mean), max(df_pos_full$rho_mean)) +
    labs(x = NULL, y = NULL) +
    # labs(x = "position score (phi)", y = "AA sensitivity scaling (rho)",
    #      title = paste("position", pos_idx)) +
    theme_classic() +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.position = "none"
    )
  cowplot::save_plot(filename = file.path(sdir_plot, "main_plot2",
                                          paste("cluster", color_idx, "_", pos_idx, ".png", sep = "")), 
                     plot = p, base_height = 3, base_width = 4, bg = "transparent")
  
}
rm(i, color_idx, pos_idx, 
   color_grey, color_blue, color_orange, color_purple,
   pos_grey, pos_blue, pos_orange, pos_purple, 
   p, df_plot_pos)


############ main plot 3: position on the structure ##########
library(bio3d)
source(file.path("workspace/newdata_0606", "structure_utils.R"))

df_pos_summarise <- df_pos_summarise %>% 
  mutate(cluster_purple = (cluster == 2) | (cluster == 3),
         cluster_blue = (cluster == 1) | (cluster == 3))
  
protein_structure <- read.pdb(file.path("workspace/met_0812/pdb", "3dkc.pdb"))

cluster1 <- map_scores_pdb(protein_structure, df_pos_summarise, "cluster_purple")
cluster2 <- map_scores_pdb(protein_structure, df_pos_summarise, "cluster_blue")
write.pdb(cluster1, file = file.path("workspace/met_0812/pdb", "3dkc_cluster_purple.pdb"))
write.pdb(cluster2, file = file.path("workspace/met_0812/pdb", "3dkc_cluster_blue.pdb"))

# set_color my_blue, [0.004, 0.451, 0.698]
# spectrum b, white my_blue, minimum = 0, maximum = 1

# set_color my_purple, [0.639, 0.412, 0.690]
# spectrum b, white my_purple, minimum = 0, maximum = 1



############ by position ##########
p1 <- ggplot(df_pos_full %>% filter(pos %in% 1059:1158), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = cond), size = 0.5) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 
p2 <- ggplot(df_pos_full %>% filter(pos %in% 1159:1258), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = cond), size = 0.5) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 
p3 <- ggplot(df_pos_full %>% filter(pos %in% 1259:1345), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = cond), size = 0.5) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 

cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot1.png"), 
                   plot = p1, base_height = 14, base_width = 22)
cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot2.png"), 
                   plot = p2, base_height = 14, base_width = 22)
cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot3.png"), 
                   plot = p3, base_height = 14, base_width = 22)
rm(p1, p2, p3)

p1 <- ggplot(df_pos_full %>% filter(pos %in% 1059:1158, genetic == "WT") %>%
               mutate(ctrl = inh == "DMSO"), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = inh, shape = ctrl), size = 1) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 
p2 <- ggplot(df_pos_full %>% filter(pos %in% 1159:1258, genetic == "WT")%>%
               mutate(ctrl = inh == "DMSO"), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = inh, shape = ctrl), size = 1) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 
p3 <- ggplot(df_pos_full %>% filter(pos %in% 1259:1345, genetic == "WT")%>%
               mutate(ctrl = inh == "DMSO"), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = inh, shape = ctrl), size = 1) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 

cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot1_WT.png"), 
                   plot = p1, base_height = 14, base_width = 22)
cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot2_WT.png"), 
                   plot = p2, base_height = 14, base_width = 22)
cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot3_WT.png"), 
                   plot = p3, base_height = 14, base_width = 22)
rm(p1, p2, p3)


p1 <- ggplot(df_pos_full %>% filter(pos %in% 1059:1158), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = genetic), size = 0.5) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 
p2 <- ggplot(df_pos_full %>% filter(pos %in% 1159:1258), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = genetic), size = 0.5) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 
p3 <- ggplot(df_pos_full %>% filter(pos %in% 1259:1345), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = genetic), size = 0.5) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 

cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot1_genetic.png"), 
                   plot = p1, base_height = 14, base_width = 22)
cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot2_genetic.png"), 
                   plot = p2, base_height = 14, base_width = 22)
cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot3_genetic.png"), 
                   plot = p3, base_height = 14, base_width = 22)
rm(p1, p2, p3)

p1 <- ggplot(df_pos_full %>% filter(pos %in% 1059:1158), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = inh), size = 0.5) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 
p2 <- ggplot(df_pos_full %>% filter(pos %in% 1159:1258), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = inh), size = 0.5) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 
p3 <- ggplot(df_pos_full %>% filter(pos %in% 1259:1345), aes(phi_mean, rho_mean)) +
  geom_point(aes(color = inh), size = 0.5) +
  theme_classic() +
  facet_wrap(vars(pos), ncol = 10) 

cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot1_inh.png"), 
                   plot = p1, base_height = 14, base_width = 22)
cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot2_inh.png"), 
                   plot = p2, base_height = 14, base_width = 22)
cowplot::save_plot(filename = file.path(sdir_plot, "pos_plot3_inh.png"), 
                   plot = p3, base_height = 14, base_width = 22)
rm(p1, p2, p3)

# ggplot(df_pos_full %>% filter(pos %in% 1219, genetic == "WT") %>%
#          mutate(ctrl = inh == "DMSO"), 
#        aes(phi_mean, rho_mean)) +
#   geom_point(aes(color = cond, shape = ctrl)) +
#   theme_classic() +
#   facet_wrap(vars(pos)) 


############ nu full ##########
df_nu_full <- data.frame()
for (cond in cond_list) {
  df_nu_full <- rbind(df_nu_full, 
                      read_tsv(file.path(dir, cond, "table", "df_nu3.tsv")) %>% 
                        mutate(cond = cond))
}
df_nu_full <- df_nu_full %>% separate(cond, into = c("gene", 'inh', "genetic"), remove = FALSE)
ggplot(df_nu_full, aes(x = as.factor(blosum_score))) +
  geom_bar(aes(y = nu_mean), stat = "identity", 
           fill = "skyblue", alpha = 0.7) +
  geom_point(aes(y = nu_mean), colour="orange", alpha = 0.9) +
  geom_errorbar(aes(ymin = nu_mean - nu_sd, ymax = nu_mean + nu_sd), 
                width=0.3, colour="orange", alpha=0.9, size=1) +
  labs(x = "blosum", y = "nu") +
  geom_text(aes(y = max(nu_mean), label = n), size = 2, vjust = -1) +
  facet_wrap(vars(cond), ncol = 6) +
  theme_classic()


