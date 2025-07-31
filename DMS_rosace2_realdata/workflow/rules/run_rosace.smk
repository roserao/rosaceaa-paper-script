rule run_rosace:
    input:
        "data/{data}.rda"
    output:
        "results/models/{data}/model{model}/rosace_eval.rds"
    log: 
        stdout = "logs/{data}/model{model}/rosace.stdout",
        stderr = "logs/{data}/model{model}/rosace.stderr"
    threads: 4
    benchmark:
        "results/models/{data}/model{model}/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R \
            --data {wildcards.data} \
            --model {wildcards.model} \
            --domain FALSE \
            > {log.stdout} 2> {log.stderr}
        """

# Rscript workflow/scripts/run_rosace.R --data metex14 --model 1 --domain FALSE

# rule all:
#     input:
#         expand("results/models/{data}/model{model}/rosace_eval.rds",
#             data = v_data, model = v_model)

rule run_rosace_domainase:
    input:
        "data/domainase_flat/{data}.rda"
    output:
        "results_domainase/models/{data}/model{model}/rosace_eval.rds"
    log: 
        stdout = "logs/domainase/{data}/model{model}/rosace.stdout",
        stderr = "logs/domainase/{data}/model{model}/rosace.stderr"
    threads: 4
    benchmark:
        "results_domainase/models/{data}/model{model}/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R \
            --data {wildcards.data} \
            --model {wildcards.model} \
            --domain TRUE \
            > {log.stdout} 2> {log.stderr}
        """

# Rscript workflow/scripts/run_rosace.R --data A0A2R8Y422_PF00240_2 --model 1 --domain TRUE    

# rule all:
#     input:
#         expand("results_domainase/models/{data}/model{model}/rosace_eval.rds",
#             data = lines, model = v_model)

rule run_rosace_met:
    input:
        "data/met2/{data}.rda"
    output:
        "results_met/models/{data}/model{model}/rosace_eval.rds"
    log: 
        stdout = "logs/met/{data}/model{model}/rosace.stdout",
        stderr = "logs/met/{data}/model{model}/rosace.stderr"
    threads: 4
    benchmark:
        "results_met/models/{data}/model{model}/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R \
            --data {wildcards.data} \
            --model {wildcards.model} \
            --met TRUE \
            > {log.stdout} 2> {log.stderr}
        """

# Rscript workflow/scripts/run_rosace.R --data met2_A458_Ex14 --model 1 --met TRUE

# rule all:
#     input:
#         expand("results_met/models/{data}/model{model}/rosace_eval.rds",
#             data = lines, model = v_model)

rule run_rosace_oct:
    input:
        "data/oct1_sm73/{data}.rda"
    output:
        "results_oct/models/{data}/model{model}/rosace_eval.rds"
    log: 
        stdout = "logs/oct/{data}/model{model}/rosace.stdout",
        stderr = "logs/oct/{data}/model{model}/rosace.stderr"
    threads: 4
    benchmark:
        "results_oct/models/{data}/model{model}/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R \
            --data {wildcards.data} \
            --model {wildcards.model} \
            --oct TRUE \
            > {log.stdout} 2> {log.stderr}
        """

# Rscript workflow/scripts/run_rosace.R --data oct1_sm73_random20_1 --model 1 --oct TRUE

# rule all:
#     input:
#         expand("results_oct/models/{data}/model{model}/rosace_eval.rds",
#             data = lines, model = v_model)

rule prelim_plot:
    input:
        "results/models/{data}/model1/rosace_eval.rds",
        "results/models/{data}/model2/rosace_eval.rds",
        "results/models/{data}/model3/rosace_eval.rds"
    output:
        "results/analysis/{data}/rosace_eval_full.rds"
    shell:
        """
        Rscript workflow/scripts/prelim_plot.R --data {wildcards.data} --domain FALSE
        """

# Rscript workflow/scripts/prelim_plot.R --data metex14 --domain FALSE
# Rscript workflow/scripts/prelim_plot.R --data oct1_sm73 --domain FALSE

# rule all:
#     input:
#         expand("results/analysis/{data}/rosace_eval_full.rds", data = v_data)

rule prelim_plot_domainase:
    input:
        "results_domainase/models/{data}/model1/rosace_eval.rds",
        "results_domainase/models/{data}/model2/rosace_eval.rds",
        "results_domainase/models/{data}/model3/rosace_eval.rds"
    output:
        "results_domainase/analysis/{data}/rosace_eval_full.rds"
    shell:
        """
        Rscript workflow/scripts/prelim_plot.R --data {wildcards.data} --domain TRUE
        """

# Rscript workflow/scripts/prelim_plot.R --data A0A2R8Y422_PF00240_2 --domain TRUE
# Rscript workflow/scripts/prelim_plot.R --data Q8IZT6_PF00612_1877 --domain TRUE

# rule all:
#     input:
#         expand("results_domainase/analysis/{data}/rosace_eval_full.rds", data = lines)

rule prelim_plot_met:
    input:
        "results_met/models/{data}/model1/rosace_eval.rds",
        "results_met/models/{data}/model2/rosace_eval.rds",
        "results_met/models/{data}/model3/rosace_eval.rds"
    output:
        "results_met/analysis/{data}/rosace_eval_full.rds"
    shell:
        """
        Rscript workflow/scripts/prelim_plot.R --data {wildcards.data} --met TRUE
        """

# Rscript workflow/scripts/prelim_plot.R --data oct1_sm73_random20_1 --met TRUE

# rule all:
#     input:
#         expand("results_met/analysis/{data}/rosace_eval_full.rds", data = lines)

rule prelim_plot_oct:
    input:
        "results_oct/models/{data}/model1/rosace_eval.rds",
        "results_oct/models/{data}/model2/rosace_eval.rds",
        "results_oct/models/{data}/model3/rosace_eval.rds"
    output:
        "results_oct/analysis/{data}/rosace_eval_full.rds"
    shell:
        """
        Rscript workflow/scripts/prelim_plot_oct.R --data {wildcards.data} 
        """

# Rscript workflow/scripts/prelim_plot_oct.R --data oct1_sm73_choice40_5

# rule all:
#     input:
#         expand("results_oct/analysis/{data}/rosace_eval_full.rds", data = lines)
