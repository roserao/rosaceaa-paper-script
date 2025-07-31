library(rosace)
library(readr)
library(stringr)
library(dplyr)
library(ggplot2)

########## directory ##########
# dir <- "data/MET2/raw_count/MET_gejf02_config_exp"
# gene <- "WT"

# dir <- "data/MET2/raw_count/MET_Ex_gejf02_config_exp"
# gene <- "Ex14"


########## read data ##########
inh_list <- c("DMSO", "A458", 
              "Crizo", "Camp", "Tepo", "Tiv", "Glu", "Savo", "NVP",
              "Cabo", "Mere", "Gle")
for (inh in inh_list) {
  create_rosace_object(inh = inh, gene = gene)
}



########## helper function ##########
func_map <- function(wt, mut) {
  if (nchar(wt) == 0) {
    return(NA)
  }
  if (wt == mut) {
    return("synonymous")
  } else if (mut == "*") {
    return("nonsense")
  } else {
    return("missense")
  }
}


create_rosace_object <- function(inh, gene) {
  key <- paste(inh, gene, sep = "_")
  
  count1 <- 
    read_tsv(file.path(dir, paste(inh, "R1_sel", sep = "_"), "main_identifiers_counts_unfiltered.tsv")) %>%
    replace(is.na(.), 0) %>%
    arrange(hgvs)
  count2 <- 
    read_tsv(file.path(dir, paste(inh, "R2_sel", sep = "_"), "main_identifiers_counts_unfiltered.tsv")) %>%
    replace(is.na(.), 0) %>%
    arrange(hgvs)
  count3 <- 
    read_tsv(file.path(dir, paste(inh, "R3_sel", sep = "_"), "main_identifiers_counts_unfiltered.tsv")) %>%
    replace(is.na(.), 0) %>%
    arrange(hgvs)
  
  var.data <- data.frame(variants = count1$hgvs)
  var.data <- var.data %>%
    mutate(tmp = substr(variants, 4, nchar(variants) - 1),
           pos = as.numeric(gsub("[[:alpha:]]", "", tmp)),
           wt = substr(tmp, 1, 1),
           tmp = substr(tmp, 2, nchar(tmp)),
           mut = gsub("[[:digit:]]", "", tmp)) %>%
    dplyr::select(-tmp)
  var.data$mut[var.data$mut == "X"] <- "*"
  var.data <- var.data %>%
    rowwise() %>%
    mutate(type = func_map(wt, mut)) %>%
    ungroup()
  var.data$ctrl <- (var.data$type == "synonymous")
  var.data$stop <- (var.data$type == "nonsense")
  
  assay1 <- CreateAssayObject(counts = as.matrix(count1[2:5]),
                              var.names = count1$hgvs,
                              key = key, rep = 1, type = "growth")
  assay2 <- CreateAssayObject(counts = as.matrix(count2[2:5]),
                              var.names = count2$hgvs,
                              key = key, rep = 2, type = "growth")
  assay3 <- CreateAssayObject(counts = as.matrix(count3[2:5]),
                              var.names = count3$hgvs,
                              key = key, rep = 3, type = "growth")
  
  rosace <- CreateRosaceObject(object = assay1, var.data = var.data)
  if (key != "A458_WT") {
    rosace <- AddAssayData(object = rosace, assay = assay2, var.data = var.data)
  }
  rosace <- AddAssayData(object = rosace, assay = assay3, var.data = var.data)
  
  rosace@var.data <- rosace@var.data %>% left_join(var.data)
  
  # preprocessing and analyze
  rosace <- NormalizeData(rosace, key = key,
                          normalization.method = "wt", 
                          wt.var.names = c("_wt"), wt.rm = TRUE)
  rosace <- IntegrateData(object = rosace, key = key)
  
  # sdir <- file.path("data", "MET2", paste("MET2", key, sep = "_"))
  # if (!dir.exists(sdir)) {
  #   dir.create(sdir)
  # }
  sdir <- file.path("data", "MET2")
  save(rosace, file = file.path(sdir, paste("met2_", key, ".rda", sep = "")))
}







