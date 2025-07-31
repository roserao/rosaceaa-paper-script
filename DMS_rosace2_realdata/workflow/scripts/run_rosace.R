#!/usr/bin/env Rscript
library('rosace')

##### parse argument #####
library('optparse')
option_list <- list(
    make_option(c("-d", "--data"), type = "character", default = NULL, 
        help = "name of protein"),
    make_option(c("-m", "--model"), type = "numeric", default = NULL, 
        help = "model mode: 1(vanilla), 2(global blosum), 3(activated blosum)"),
    make_option(c("-f", "--domain"), type = "logical", default = FALSE, 
        help = "whether to process domainase data"),
    make_option(c("-g", "--met"), type = "logical", default = FALSE, 
        help = "whether to process met2 data"),
    make_option(c("-c", "--oct"), type = "logical", default = FALSE, 
        help = "whether to process oct1 data")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data <- opt$data
idx_model <- opt$model
flag_domain <- opt$domain
flag_met <- opt$met
flag_oct <- opt$oct

# data <- "slco1b1"
# idx_model <- 1

##### load data #####
##### key ##### 
if (flag_domain) {
  load(file.path('data', 'domainase_flat', paste(data, ".rda", sep = "")))
  key <- data
} else if (flag_met) {
  load(file.path('data', "met2", paste(data, ".rda", sep = "")))
  key <- substr(data, 6, nchar(data))
} else if (flag_oct) {
  load(file.path('data', "oct1_sm73", paste(data, ".rda", sep = "")))
  key <- "1SM73"
} else {
  load(file.path('data', paste(data, ".rda", sep = "")))
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

##### save directory #####
if (flag_domain) {
  sdir <- file.path('results_domainase/models', data, paste("model", idx_model, sep = ""))
} else if (flag_met) {
  sdir <- file.path('results_met/models', data, paste("model", idx_model, sep = ""))
} else if (flag_oct) {
  sdir <- file.path('results_oct/models', data, paste("model", idx_model, sep = ""))
} else {
  sdir <- file.path('results/models', data, paste("model", idx_model, sep = ""))
}
if (!dir.exists(sdir)) {
    dir.create(sdir, recursive = TRUE)
} 

##### run Rosace #####
if (idx_model == 1) {
    # Model 1
    rosace <- RunRosace(object = rosace,
                        name = key,
                        type = "AssaySet",
                        savedir = sdir,
                        pos.col = "pos",
                        ctrl.col = "ctrl",
                        ctrl.name = TRUE,
                        stop.col = "stop",
                        stop.name = TRUE,
                        install = FALSE)
} else if (idx_model == 2) {
    # Model 2
    rosace <- RunRosace(object = rosace,
                        name = key,
                        type = "AssaySet",
                        savedir = sdir,
                        pos.col = "pos",
                        ctrl.col = "ctrl",
                        ctrl.name = TRUE,
                        stop.col = "stop",
                        stop.name = TRUE,
                        wt.col = "wt",
                        mut.col = "mut",
                        aa.code = "single",
                        pos.act = FALSE,
                        install = FALSE)

} else if (idx_model == 3) {
    # Model 3
    rosace <- RunRosace(object = rosace,
                        name = key,
                        type = "AssaySet",
                        savedir = sdir,
                        pos.col = "pos",
                        ctrl.col = "ctrl",
                        ctrl.name = TRUE,
                        stop.col = "stop",
                        stop.name = TRUE,
                        wt.col = "wt",
                        mut.col = "mut",
                        aa.code = "single",
                        pos.act = TRUE,
                        install = FALSE)
} else {
    stop("Invalid rosace model mode: 1(vanilla), 2(global blosum), 3(activated blosum)")
}

##### save Rosace object #####
saveRDS(rosace, file = file.path(sdir, "rosace_eval.rds"))
