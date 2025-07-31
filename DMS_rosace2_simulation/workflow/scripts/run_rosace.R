#!/usr/bin/env Rscript
library('rosace')

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
        help = "simulation index"),
    make_option(c("-l", "--model"), type = "numeric", default = NULL, 
        help = "model mode: 1(vanilla), 2(global blosum), 3(activated blosum)")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data_mode <- opt$data
sim_mode <- opt$mode
n_rep <- opt$rep
n_round <- opt$round
idx_sim <- opt$sim
idx_model <- opt$model

# data_mode <- 1
# sim_mode <- 1
# idx_sim <- 1
# idx_rep <- 3
# idx_round <- 4
# idx_model <- 1

##### load data #####
dir <-  
    file.path('results/counts', 
              paste('data', data_mode, sep = ""),
              paste("growth_rep", n_rep, "_rd", n_round, "_", sim_mode, sep = ""),
              paste("sim", idx_sim, sep = ""))
ddir <- file.path(dir, "rosace")
rosace <- readRDS(file.path(ddir, "rosace.rds"))
key <- "simulation"

##### save directory #####
sdir <-   
    file.path('results/models', 
              paste('data', data_mode, sep = ""),
              paste("growth_rep", n_rep, "_rd", n_round, "_", sim_mode, sep = ""),
              paste("sim", idx_sim, sep = ""),
              paste("model", idx_model, sep = ""))
if (!dir.exists(sdir)) {
    dir.create(sdir, recursive = TRUE)
} 

##### preprocessing #####
# rosace <- FilterData(rosace, key = key, na.rmax = 0.5)
rosace <- ImputeData(rosace, key = key, impute.method = "knn", na.rmax = 0.5)
rosace <- NormalizeData(rosace, 
                        key = key,
                        normalization.method = "wt", 
                        wt.var.names = rosace@var.data[[1]][rosace@var.data$ctrl], 
                        wt.rm = FALSE)
rosace <- IntegrateData(object = rosace, key = key)

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

