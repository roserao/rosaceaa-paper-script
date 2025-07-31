rule generate_counts:
    """
    'data', 'mode', 'rep', and 'rd' are wildcards.
    'sim' should be range(n_sim).
    """
    input: 
        "data/slco1b1.rda",
        "data/oct1_sm73.rda",
        "data/metex14.rda",
        "data/slco1b1_ros.rda"
    output: 
        expand("results/counts/data{{data}}/growth_rep{{rep}}_rd{{rd}}_{{mode}}/sim{sim}/rosace/rosace.rds",
            sim = range(1, 1 + nsim))
    shell: 
        """
        Rscript workflow/scripts/generate_counts.R \
            --mode {wildcards.mode} \
            --data {wildcards.data} \
            --rep {wildcards.rep} \
            --round {wildcards.rd} \
            --sim {nsim}
        """

# rule all:
#     input:
#         expand("results/counts/data{data}/growth_rep{rep}_rd{rd}_{mode}/sim{sim}/rosace/rosace.rds",
#             sim = range(1, 1 + nsim), data = v_data, rep = v_nrep, rd = v_nround, mode = v_mode)


rule run_rosace:
    input:
        "results/counts/data{data}/growth_rep{rep}_rd{rd}_{mode}/sim{sim}/rosace/rosace.rds"
    output:
        "results/models/data{data}/growth_rep{rep}_rd{rd}_{mode}/sim{sim}/model{model}/rosace_eval.rds"
    log: 
        stdout = "logs/data{data}/growth_rep{rep}_rd{rd}_{mode}/sim{sim}/model{model}/rosace.stdout",
        stderr = "logs/data{data}/growth_rep{rep}_rd{rd}_{mode}/sim{sim}/model{model}/rosace.stderr"
    threads: 4
    benchmark:
        "results/models/data{data}/growth_rep{rep}_rd{rd}_{mode}/sim{sim}/model{model}/benchmark.txt"
    shell:
        """
        Rscript workflow/scripts/run_rosace.R \
            --mode {wildcards.mode} \
            --data {wildcards.data} \
            --rep {wildcards.rep} \
            --round {wildcards.rd} \
            --sim {wildcards.sim} \
            --model {wildcards.model} \
            > {log.stdout} 2> {log.stderr}
        """

# Rscript workflow/scripts/run_rosace.R --mode 1 --data 1 --rep 3 --round 3 --sim 1 --model 1

# rule all:
#     input:
#         expand("results/models/data{data}/growth_rep{rep}_rd{rd}_{mode}/sim{sim}/model{model}/rosace_eval.rds",
#             sim = range(1, 1 + nsim), data = v_data, rep = v_nrep, rd = v_nround, mode = v_mode, model = v_model)

rule analyze_rosace:
    input:
        expand("results/models/data{{data}}/growth_rep{{rep}}_rd{{rd}}_{{mode}}/sim{{sim}}/model{model}/rosace_eval.rds",
            model = range(1, 4))
    output:
        "results/analysis/data{data}/growth_rep{rep}_rd{rd}_{mode}/sim{sim}/rosace_eval_full.rds"
    shell:
        """
        Rscript workflow/scripts/analyze_rosace.R \
            --mode {wildcards.mode} \
            --data {wildcards.data} \
            --rep {wildcards.rep} \
            --round {wildcards.rd} \
            --sim {wildcards.sim} 
        """

# Rscript workflow/scripts/analyze_rosace.R --mode 1 --data 1 --rep 3 --round 3 --sim 1 