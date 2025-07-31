library(tidyverse)
library(cowplot)

dir <- file.path("workspace", "domainase_0614")
ddir <- file.path(dir, "results_domainase", "analysis")
sdir <- file.path(dir, "results_domainase", "summary")

########## extract secondary structure ########## 

##### read secondary structure #####
df <- read.table("data/domainase/analysis_files/domainome_structural_annotation.txt")
df <- df %>% select(dom_ID, pos, pos_in_uniprot, variant_ID, 
                    WT, wt_aa, mut_aa, 
                    secondary_structure_code, secondary_structure)
df <- df %>% arrange(dom_ID)
df <- df %>% 
  select(dom_ID, pos, pos_in_uniprot, secondary_structure_code, secondary_structure) %>% 
  unique()
write_tsv(df, file = file.path(dir, "analysis_data", "secondary_structure_annotation.tsv"))

table(df$secondary_structure)
# 310Helix AlphaHelix     Bridge       Coil     Strand       Turn 
# 1432      19726        751      11394       9811      11674 
df_domain <- df %>% group_by(dom_ID) %>% summarise(length = max(pos))
ggplot(df_domain, aes(length)) +
  geom_histogram(color = "grey") +
  theme_classic()

##### combine with phi and rho score #####
length(unique(df_domain$dom_ID))
length(list.files(ddir))

df_rosace_blosum <- data.frame()
for (dom in list.files(ddir)) {
  load(file.path(ddir, dom, "df_pos_full.rda"))
  df_rosace_blosum <- rbind(
    df_rosace_blosum,
    df_pos3 %>% mutate(dom_ID = dom)
  )
}
rm(dom, df_pos1, df_pos2, df_pos3)

df_rosace_blosum <- df_rosace_blosum %>% left_join(df, by = join_by(pos, dom_ID))
write_tsv(df_rosace_blosum, file.path(dir, "analysis_data", "secondary_structure_annotation_rosace_blosum.tsv"))

##### Need a filter for significant nu #####
df_nu3 <- read_tsv(file = file.path(dir, "results_domainase", "summary", "df_nu3.tsv"))

df_nu3 <- df_nu3 %>%
  rowwise() %>%
  mutate(lfsr.neg = stats::pnorm(0, mean = nu_mean, sd = nu_sd, lower.tail = FALSE),
         lfsr.pos = stats::pnorm(0, mean = nu_mean, sd = nu_sd, lower.tail = TRUE),
         lfsr = min(.data$lfsr.neg, .data$lfsr.pos)) %>%
  ungroup()
df_nu3 <- df_nu3 %>% mutate(nu.sig = (lfsr <= 0.1))

df_nu3_domain <- df_nu3 %>% group_by(dom_ID) %>% summarise(nu.num.sig = sum(nu.sig))
df_nu3_domain <- df_nu3_domain %>% mutate(nu.sig = (nu.num.sig >= 10))

df_rosace_blosum <- df_rosace_blosum %>% left_join(df_nu3_domain)
# write_tsv(df_rosace_blosum, file.path(dir, "analysis_data", "secondary_structure_annotation_rosace_blosum_nufilter.tsv"))

##### domain meta-data #####
df_nu3_domain <- df_nu3_domain %>% left_join(df_domain)
rm(df_domain)
df_nu3_domain <- df_nu3_domain %>% separate(dom_ID, c("gene", "PFAM_ID", "fpos"), remove = FALSE, sep = "_")  
df_scop <- read_tsv("data/domainase/analysis_files/PFAM_ID_to_SCOP_class.tsv")
df_nu3_domain <- df_nu3_domain %>% left_join(df_scop)
write_tsv(df_nu3_domain, file.path(dir, "analysis_data", "domID_annotation.tsv"))

df_rosace_blosum <- df_rosace_blosum %>% left_join(df_nu3_domain %>% select(dom_ID, scop_class, zinc_finger = `zinc finger`))
write_tsv(df_rosace_blosum, file.path(dir, "analysis_data", "secondary_structure_annotation_rosace_blosum_nufilter.tsv"))


df_rosace_blosum <- read_tsv(file.path(dir, "analysis_data", "secondary_structure_annotation_rosace_blosum_nufilter.tsv"))
df_nu3_domain <- read_tsv(file.path(dir, "analysis_data", "domID_annotation.tsv"))
range(df_rosace_blosum$phi_mean)

sse_mtx <- read_tsv(file.path("workspace/newdata_0606", "secondary_structure", "sse_ABCG2_mtx.tsv"))
range(sse_mtx$phi_mean)


##### variant dms data #####
df_rosace_blosum <- data.frame()
for (dom in list.files(ddir)) {
  load(file.path(ddir, dom, "df_var_full.rda"))
  df_rosace_blosum <- rbind(
    df_rosace_blosum,
    df_var3 %>% mutate(dom_ID = dom)
  )
}
rm(dom, df_var1, df_var2, df_var3)

df_rosace_dms <- df_rosace_blosum %>% group_by(dom_ID, pos, pos_in_uniprot) %>%
  pivot_wider(id_cols = c("dom_ID", "pos", "pos_in_uniprot"),
               names_from = "mut", values_from = "mean")
df_rosace_blosum_pos <- read_tsv(file.path(dir, "analysis_data", "secondary_structure_annotation_rosace_blosum.tsv"))
df_rosace_blosum_pos <- df_rosace_blosum_pos %>% left_join(df_rosace_dms)
write_tsv(df_rosace_blosum_pos, file.path(dir, "analysis_data", "secondary_structure_annotation_rosace_blosum_variant.tsv"))

