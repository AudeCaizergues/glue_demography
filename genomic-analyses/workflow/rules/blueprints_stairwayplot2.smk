# Rule for generating the blueprint files to run stairwayplot2 for estimations of Ne

#######################
### BLUEPRINT FILES ###
#######################

rule create_blueprints_files_stairwayplot:
    """
    Create blueprint files from SFS, for stairwayplot2.
    """
    input:
        sfs = rules.angsd_estimate_sfs_byCity_byHabitat.output,
        bams = rules.create_bam_list_byCity_byHabitat_withoutRelated.output,
        logfile = rules.angsd_estimate_sfs_byCity_byHabitat.log,
    output:
        '{0}/{{city}}_{{habitat}}_{{site}}.blueprint'.format(STAIRWAYPLOT_DIR)
    log: LOG_DIR + '/blueprint_stairwayplot/{city}_{habitat}_{site}.blueprint.log'
    shell:
        """
          city_name={wildcards.city}
          habitat_name={wildcards.habitat}
          echo "popid: ""$city_name""_""$habitat_name" > {output}
          nseq=$(( $(wc -l {input.bams} | cut -d " " -f 1) * 2 ))
          echo "nseq: ""$nseq" >> {output}
          echo "L: ""$(grep "Will run optimization on nSites:" {input.logfile} | cut -d " " -f 7)" >> {output}
          echo "wether_folded: true" >> {output}
          echo "SFS: ""$(head -n1 {input.sfs} | sed 's/ /	/g' |cut -d $'\t' -f 2-$((($nseq+2)/2))) " >> {output}
          echo "pct_training: 0.67" >> {output}
          echo "nrand: 20 40 60 80" >> {output}
          echo "project_dir: /scratch/projects/trifolium/glue/demography/glue_demography/results/stairwayplot/byCity_byHabitat/""$city_name""/""$habitat_name""/" >> {output}
          echo "stairway_plot_dir: /scratch/projects/trifolium/glue/demography/glue_demography/genomic-analyses/workflow/scripts/stairwayplot2/stairway_plot_v2.1.1/stairway_plot_es/" >> {output}
          echo "ninput: 600" >> {output}
          echo "random_seed: 42" >> {output}
          echo "mu: 1.2e-8" >> {output}
          echo "year_per_generation: 1" >> {output}
          echo "plot_title: ""$city_name ""$habitat_name" >> {output}
          mkdir -p /scratch/projects/trifolium/glue/demography/glue_demography/results/stairwayplot/byCity_byHabitat/{wildcards.city}/{wildcards.habitat}/
        """

rule run_blueprints_files_stairwayplot:
    """
    Run first stairwayplot2 step.
    """
    input:
        blueprints = rules.create_blueprints_files_stairwayplot.output,
    output:
        '{0}/{{city}}_{{habitat}}_{{site}}.blueprint'.format(STAIRWAYPLOT_DIR)
    log: LOG_DIR + '/blueprint_stairwayplot/{city}_{habitat}_{site}.blueprint.log'
    shell:
        """
          java -cp stairway_plot_es Stairbuilder 
        """

##############
#### POST ####
##############

rule create_blueprints_done:
    input:
        expand(rules.create_blueprints_files_stairwayplot.output, habitat= HABITATS, city=CITIES, site='4fold', maf='0.05')
    output:
        '{0}/blueprints.done'.format(STAIRWAYPLOT_DIR)
    shell:
        """
        touch {output}
        """