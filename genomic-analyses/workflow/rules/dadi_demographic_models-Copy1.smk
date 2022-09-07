# Rules to perform demographic modelling using dadi


rule run_dadi:
    input:
        sfs = rules.format_dadi_sfs.output
    output:
        logf = '{0}/{{city}}_pop0_pop1_{{rep}}.{{model}}.log.txt'.format(DADI_DIR),
        optimf = '{0}/{{city}}_pop0_pop1_{{rep}}.{{model}}.optimized.txt'.format(DADI_DIR)
    log: LOG_DIR + '/run_dadi/{city}_{model}_{rep}.log'
    conda: '../envs/dadi.yaml'
    params:
        prefix = '{0}/{{model}}/'.format(DADI_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    script:
        "../scripts/python/dadi_pipeline/dadi_Run_2D_Set.py"

rule summarize_dadi_output:
    input:
        expand(rules.run_dadi.output, city=CITIES, model=DADI_MODELS, rep = ['1', '2', '3', '4', '5'])
    output:
        short = '{0}/{{city}}/Results_Summary_Short.txt'.format(DADI_DIR),
        exten = '{0}/{{city}}/Results_Summary_Extended.txt'.format(DADI_DIR)
    params:
        path = '{0}/'.format(DADI_DIR)
    shell:
        """
        python3 scripts/python/dadi_pipeline/Summarize_Outputs.py {params.path}
        """

rule dadi_done:
    input:
        rules.summarize_dadi_output.output
    output:
        '{0}/dadi.done'.format(DADI_DIR)
    shell:
        """
        touch {output}
        """
