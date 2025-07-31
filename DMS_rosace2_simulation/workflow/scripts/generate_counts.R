#!/usr/bin/env Rscript
library('rosace')

##### parse argument #####
library('optparse')
option_list <- list(
    make_option(c("-m", "--mode"), type = "numeric", default = NULL, 
        help = "simualation mode: 1, 2, 3, 4, 5"),
    make_option(c("-d", "--data"), type = "numeric", default = NULL, 
        help = "data code: 1(slco1b1), 2(oct1+sm73), 3(metex14), 4(slco1b1+ros)"),
    make_option(c("-s", "--sim"), type = "numeric", default = NULL, 
        help = "number of simulations"),
    make_option(c("-r", "--rep"), type = "numeric", default = NULL,
        help = "number of replicates"),
    make_option(c("-t", "--round"), type = "numeric", default = NULL, 
        help = "number of rounds")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
data_mode <- opt$data
sim_mode <- opt$mode
n_sim <- opt$sim
n_rep <- opt$rep
n_round <- opt$round

# data_mode <- 1
# sim_mode <- 1
# n_sim <- 5
# n_rep <- 3
# n_round <- 4

##### generate config from rosace #####
# please test parameter
if (data_mode == 1) {
    load("data/slco1b1.rda")
    cfg <- 
        CreateConfig(rosette,
                    n.sim = n_sim, mode.sim = sim_mode,
                    n.rep = n_rep, n.round = n_round, 
                    type.sim = "growth", save.sim = "results/counts/data1/",
                    lib.shrink = 1, seq.shrink = 1, 
                    wt.effect = 0.5, pop.size = 100, seq.depth = 200)
} else if (data_mode == 2) {
    load("data/oct1_sm73.rda")
    cfg <- 
        CreateConfig(rosette,
                    n.sim = n_sim, mode.sim = sim_mode,
                    n.rep = n_rep, n.round = n_round, 
                    type.sim = "growth", save.sim = "results/counts/data2/",
                    lib.shrink = 1, seq.shrink = 1, 
                    wt.effect = 1, pop.size = 500, seq.depth = 100)
} else if (data_mode == 3) {
    load("data/metex14.rda")
    cfg <- 
        CreateConfig(rosette,
                    n.sim = n_sim, mode.sim = sim_mode,
                    n.rep = n_rep, n.round = n_round, 
                    type.sim = "growth", save.sim = "results/counts/data3/",
                    lib.shrink = 1, seq.shrink = 1, var.shrink = 1,
                    wt.effect = -1, pop.size = 200, seq.depth = 300)
} else if (data_mode == 4) {
    load("data/slco1b1_ros.rda")
    cfg <- 
        CreateConfig(rosette,
                    n.sim = n_sim, mode.sim = sim_mode,
                    n.rep = n_rep, n.round = n_round, 
                    type.sim = "growth", save.sim = "results/counts/data4/",
                    lib.shrink = 1, seq.shrink = 1, 
                    wt.effect = 0.5, pop.size = 100, seq.depth = 200)
} else {
    stop("Invalid data code. Please enter 1-4: 1(slco1b1), 2(oct1+sm73), 3(metex14), 4(slco1b1+ros)")
}

##### run Rosette simulation #####
runRosette(config = cfg, save.tsv = TRUE, save.rosace = TRUE, save.enrich2 = FALSE)
