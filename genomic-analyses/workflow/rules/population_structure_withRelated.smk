# Population structure analyses. 

#######################
### PCA & ADMIXTURE ###
#######################

# PCA & admixture analysis using LD-pruned 4fold SNPs from above (MAF > 0.05)

####################
#### GLOBAL PCA ####
####################

rule pcangsd_allSamples_withRelated:
    """
    Perform PCA using genome-wide 4fold dengenerate sites using all samples from all cities.
    """
    input:
        rules.concat_angsd_gl_pruned.output
    output:
        '{0}/pcangsd/allSamples/pcangsd_withRelated/allSamples_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
    log: LOG_DIR + '/pcangsd_allSamples/pcangsd_withRelated/allSamples_{site}_maf{maf}_pcangsd.log'
    container: 'library://james-s-santangelo/pcangsd/pcangsd:0.99'
    threads: 12
    params:
        out = '{0}/pcangsd/allSamples/pcangsd_withRelated/allSamples_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 12000,
        time = '03:00:00'
    shell:
        """
        python3 /opt/pcangsd-v.0.99/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -minMaf {wildcards.maf} \
            -threads {threads} \
            -iter 5000 \
            > {log}
        """

#####################
#### PCA BY CITY ####
#####################


rule pcangsd_byCity_withRelated:
    """
    Perform PCA by city using genome-wide 4fold dengenerate sites and estimate admixture proportions
    """
    input:
        rules.angsd_gl_byCity.output.gls
    output:
        cov = '{0}/pcangsd/by_city/pcangsd_withRelated/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
        Q = '{0}/pcangsd/by_city/pcangsd_withRelated/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd.admix.Q.npy'.format(POP_STRUC_DIR)
    log: LOG_DIR + '/pcangsd_byCity/pcangsd_withRelated/{city}_{site}_maf{maf}_pcangsd.log'
    container: 'library://james-s-santangelo/pcangsd/pcangsd:0.99'
    threads: 12
    params:
        out = '{0}/pcangsd/by_city/pcangsd_withRelated/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '02:00:00'
    shell:
        """
        python3 /opt/pcangsd-v.0.99/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -minMaf {wildcards.maf} \
            -threads {threads} \
            -admix \
            -admix_seed 42 \
            -iter 5000 \
            > {log}
        """

        
##############
#### POST ####
##############

rule pop_structure_withRelated_done:
    """
    Generate empty flag file signaling successful completion of PCAngsd
    """
    input:
        expand(rules.pcangsd_allSamples_withRelated.output, site = '4fold', maf = ['0.05']),
        expand(rules.pcangsd_byCity_withRelated.output, site = '4fold', maf = '0.05', city = CITIES)
    output:
        '{0}/population_structure_withRelated.done'.format(POP_STRUC_DIR)
    shell:
        """
        touch {output}
        """
