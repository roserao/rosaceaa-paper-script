library(rosace)
library(tidyverse)

df <- read_table("data/domainase/analysis_files/mutated_domainome_merged_filtered_all_VEPs.txt")
df_count <- data.frame(table(df$dom_ID))

df_count2 <- data.frame(table(df$dom_ID, df$library)) %>% filter(Freq != 0)
table(df_count2$Var2)
df_count2 <- df_count2 %>% arrange(Var1)
write.table(df_count2$Var1, 
            file.path(sdir, "domain_list.txt"),
            quote = FALSE, col.names = FALSE, row.names = FALSE)


ddir1 <- "data/domainase/dimsum_scripts_and_inputfiles"
ddir2 <- "data/domainase/dimsum_output"
sdir <- "workspace/domainase_0614/data"

domain_remove <- c()
domain_totalnorm <- c()
code_list <- c("A1", "B3", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
for (code in code_list) {
  
  ssdir <- file.path(sdir, code)
  if(!dir.exists(ssdir)) {
    dir.create(ssdir)
  }
  
  df1 <- read_table(file.path(ddir1, paste(code, "all_aa_variants.txt", sep = "_")))
  if (code  == "A1") {
    load(file.path(ddir2, paste(code, "BGI_Q30_fitness_replicates.RData", sep = "_")))
  } else {
    load(file.path(ddir2, paste(code, "BGI_Q20_fitness_replicates.RData", sep = "_")))
  }
  rm(doubles, singles, wildtype, synonymous)

  df1 <- df1 %>% 
    left_join(all_variants %>% select(aa_seq, starts_with("count"), fitness, sigma),
              by = c("aa_seq" = "aa_seq"))
  rm(all_variants)
  
  df1 <- df1 %>% rowwise() %>%
    filter(!all(is.na(count_e1_s0), is.na(count_e1_s1),
                is.na(count_e2_s0), is.na(count_e2_s1),
                is.na(count_e3_s0), is.na(count_e3_s1))) %>%
    ungroup()
  
  domain_list <- as.character(df_count2[df_count2$Var2 == code, ]$Var1)
  
  for (domain in domain_list) {
    df_sub <- df1 %>% filter(dom_ID == domain)
    type <- "growth"
    
    if (domain == "A1X283_PF00018_155") {
      df_sub <- df_sub %>% filter(wt_seq == "EQYVVVANYQKQESSEISLSVGQVVDIIEKNESGWWFVSTAEEQGWVPATCLEGQDGV")
    }
    
    if (nrow(df_sub) < 10) {
      warning(paste("not enough variant for protein domain", domain, sep = " "))
      domain_remove <- c(domain_remove, domain)
      next
    }
    
    df_sub <- df_sub %>%
      unite("variants", wt_aa, pos, mut_aa, remove = FALSE) %>%
      relocate(variants)
    
    assay1 <- CreateAssayObject(counts = as.matrix(df_sub[c(11, 14)]),
                                var.names = df_sub$variants,
                                key = domain, rep = 1, type = type)
    assay2 <- CreateAssayObject(counts = as.matrix(df_sub[c(12, 15)]),
                                var.names = df_sub$variants,
                                key = domain, rep = 2, type = type)
    assay3 <- CreateAssayObject(counts = as.matrix(df_sub[c(13, 16)]),
                                var.names = df_sub$variants,
                                key = domain, rep = 3, type = type)
    
    rosace <- CreateRosaceObject(object = assay1)
    rosace <- AddAssayData(object = rosace, assay = assay2)
    rosace <- AddAssayData(object = rosace, assay = assay3)
    GetAssayName(rosace)
    
    rosace <- FilterData(rosace, key = domain, na.rmax = 0.5, min.count = 20)
    rosace <- ImputeData(rosace, key = domain, impute.method = "zero")
    
    if (any(df_sub$WT)) {
      rosace <- NormalizeData(rosace, key = domain,
                              normalization.method = "wt", 
                              wt.var.names = df_sub$variants[df_sub$WT], wt.rm = FALSE)
    } else {
      rosace <- NormalizeData(rosace, key = domain, normalization.method = "total")
      domain_totalnorm <- c(domain_totalnorm, domain)
    }
    
    rosace <- IntegrateData(object = rosace, key = domain)
    GetAssaySetName(rosace)
    
    rosace@var.data <- df_sub[c(1:5, 7:10)]
    rosace@var.data$ctrl <- df_sub$WT
    rosace@var.data$stop <- (df_sub$mut_aa == "*")
    
    rosace@var.data <- rosace@var.data %>% 
      rename(pos = pos, wt = wt_aa, mut = mut_aa)
    
    save(rosace, file = file.path(ssdir, paste(domain, ".rda", sep = "")))
  }
  
  write.table(domain_remove, 
              file.path(sdir, paste(code, "remove.txt", sep = "_")),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  write.table(domain_totalnorm, 
              file.path(sdir, paste(code, "totalnorm.txt", sep = "_")),
              quote = FALSE, col.names = FALSE, row.names = FALSE)
  domain_remove <- c()
  domain_totalnorm <- c()
}
