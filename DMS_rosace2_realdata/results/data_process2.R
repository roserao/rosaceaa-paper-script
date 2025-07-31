library(rosace)
library(readr)
library(tidyr)
library(dplyr)


########## new data list ##########
dir <- "workspace/review2data_0602_25/analysis"
cond_list <- list.dirs(dir, full.names = FALSE, recursive = FALSE)

##### load variant data frame and compute summary statistics #####
df_decomp <- data.frame()
for (cond in cond_list) {
  load(file.path(dir, cond, "df_var_full.rda"))
  df_decomp_key <- data.frame(
    dom_ID = cond,
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
rm(cond, df_decomp_key)

df_decomp %>% filter(model == 3)
