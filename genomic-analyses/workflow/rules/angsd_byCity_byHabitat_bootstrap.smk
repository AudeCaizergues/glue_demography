# Rules for estimating bootstraped SFS (1D and 2D), summary stats, and GLs for urban and rural habitats within cities

###############################
#### SFS AND SUMMARY STATS ####
###############################

rule angsd_estimate_joint_sfs_byCity_bootstrap:
    """
    Estimated folded, two-dimensional urban-rural SFS for each city using realSFS. Uses 4fold sites.
    """
    input:
        saf = ancient(get_habitat_saf_files_byCity),
        sites = rules.convert_sites_for_angsd.output,
        idx = rules.angsd_index_degenerate_sites.output,
    output:
        '{0}/sfs/by_city/bootstrap/{{city}}/{{city}}_{{site}}_r_u.2dsfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_2dsfs_byCity/bootstrap/{city}_{site}.2dsfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '03:00:00'
    shell:
        """
        realSFS {input.saf} \
            -sites {input.sites} \
            -maxIter 2000 \
            -seed 42 \
            -fold 1 \
            -bootstrap 100 \
            -P {threads} > {output} 2> {log}
        """

rule angsd_estimate_sfs_byCity_byHabitat_bootstrap:
    """
    Estimate folded SFS separately for each habitat in each city (i.e., 1D SFS) using realSFS woth 100 bootstrap. 
    """
    input:
        saf = ancient(rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx),
        sites = rules.convert_sites_for_angsd.output,
        idx = rules.angsd_index_degenerate_sites.output,
    output:
        '{0}/sfs/by_city/bootstrap/{{city}}/{{city}}_{{habitat}}_{{site}}.sfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_sfs_byCity_byHabitat/bootstrap/{city}_{habitat}_{site}_sfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 5000,
        time = '03:00:00'
    shell:
        """
        realSFS {input.saf} \
            -sites {input.sites} \
            -P {threads} \
            -fold 1 \
            -maxIter 2000 \
            -bootstrap 100 \
            -seed 42 > {output} 2> {log}
        """

rule split_bootstrapped_sfs_byCity:
    """
    Estimated folded, two-dimensional urban-rural SFS for each city using realSFS. Uses 4fold sites.
    """
    input:
        ancient(rules.angsd_estimate_joint_sfs_byCity_bootstrap.output)
    output:
        '{0}/sfs/by_city/bootstrap/{{city}}/{{city}}_{{site}}_r_u_split.2dsfs{{number}}'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_2dsfs_byCity/bootstrap/split_bootstrapped_{city}_{site}_{{number}}.2dsfs.log'
    params : 
        out = '{0}/sfs/by_city/bootstrap/{{city}}/{{city}}_{{site}}_r_u_split.2dsfs'.format(ANGSD_DIR)
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '03:00:00'
    shell:
        """
        split -d -l 1 {input} {params.out} 2> {log}
        """
        
rule split_bootstrapped_sfs_byCity_byHabitat:
    """
    Estimated folded, two-dimensional urban-rural SFS for each city using realSFS. Uses 4fold sites.
    """
    input:
        rules.angsd_estimate_sfs_byCity_byHabitat_bootstrap.output
    output:
        '{0}/sfs/by_city/bootstrap/{{city}}/{{city}}_{{habitat}}_{{site}}_split.sfs{{number}}'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_2dsfs_byCity/bootstrap/split_bootstrapped_{city}_{habitat}_{site}_{{number}}.sfs.log'
    params : 
        out = '{0}/sfs/by_city/bootstrap/{{city}}/{{city}}_{{habitat}}_{{site}}_split.sfs'.format(ANGSD_DIR)
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '03:00:00'
    shell:
        """
        split -d -l 1 {input} {params.out} 2> {log}
        """

#####################################
#### BOOTSTRAPPED FST AND THETAS ####
#####################################

rule angsd_bootsrapped_fst_index:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Hudson's Fst estimator.
    Uses bootstrapped urban and rural samples.
    """
    input: 
        saf_idx = get_habitat_saf_files_byCity_bootstrap,
        joint_sfs = rules.split_bootstrapped_sfs_byCity.output,
        sites = rules.convert_sites_for_angsd.output
    output:
        fst = '{0}/summary_stats/hudson_fst/bootstrap/{{city}}/{{city}}_{{site}}_{{number}}_r_u.fst.gz'.format(ANGSD_DIR),
        idx = '{0}/summary_stats/hudson_fst/bootstrap/{{city}}/{{city}}_{{site}}_{{number}}_r_u.fst.idx'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_bootstrapped_fst_index/{city}_{site}_{number}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    params:
        fstout = '{0}/summary_stats/hudson_fst/{{city}}/bootstrap/{{city}}_{{site}}_{{number}}_r_u'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf_idx} \
            -sites {input.sites} \
            -sfs {input.joint_sfs} \
            -fold 1 \
            -P {threads} \
            -whichFst 1 \
            -fstout {params.fstout} 2> {log}
        """

rule angsd_bootstrapped_fst_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. Uses permuted 
    urban and rural samples.
    """
    input:
        rules.angsd_bootsrapped_fst_index.output.idx
    output:
        '{0}/summary_stats/hudson_fst/{{city}}/bootstrap/{{city}}_{{site}}_{{number}}_r_u_readable.fst'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_bootstrapped_fst_readable/{city}_{site}_{number}_readable.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

rule angsd_estimate_bootstrapped_thetas_byCity_byHabitat:
    """
    Generate per-site thetas in each habitat for each city from 1DSFS. Uses permuted urban and rural samples.
    """
    input:
        saf_idx = ancient(rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx),
        sfs = rules.split_bootstrapped_sfs_byCity_byHabitat.output,
        sites = rules.convert_sites_for_angsd.output
    output:
        idx = '{0}/summary_stats/thetas/by_city/{{city}}/bootstrap/{{city}}_{{habitat}}_{{site}}_{{number}}.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/by_city/{{city}}/bootstrap/{{city}}_{{habitat}}_{{site}}_{{number}}.thetas.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_bootstrapped_thetas_byCity_byHabitat/{city}_{habitat}_{site}_{number}_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 4
    params:
        out = '{0}/summary_stats/thetas/by_city/{{city}}/bootstrap/{{city}}_{{habitat}}_{{site}}_{{number}}'.format(ANGSD_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        realSFS saf2theta {input.saf_idx} \
            -sites {input.sites} \
            -P {threads} \
            -fold 1 \
            -sfs {input.sfs} \
            -outname {params.out} 2> {log}
        """

rule angsd_bootstrapped_diversity_neutrality_stats_byCity_byHabitat:
    """
    Estimate pi, Waterson's theta, Tajima's D, etc. in each habitat in each city. Uses permuted urban and rural samples
    """
    input:
        rules.angsd_estimate_bootstrapped_thetas_byCity_byHabitat.output.idx
    output:
       '{0}/summary_stats/thetas/by_city/{{city}}/bootstrap/{{city}}_{{habitat}}_{{site}}_{{number}}.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_bootstrapped_diversity_neutrality_stats_byCity_byHabitat/{city}_{habitat}_{site}_{number}_diversity_neutrality.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """

##############
#### POST ####
##############

rule angsd_byCity_byHabitat_bootstrap_done:
    """
    Generate empty flag file signalling successful completion of SFS, summary stat and GL estimation 
    for habitats within cities
    """
    input:
        expand(rules.angsd_estimate_joint_sfs_byCity_bootstrap.output, city=CITIES, site=['4fold']),
        expand(rules.split_bootstrapped_sfs_byCity.output, city=CITIES, number=NUMBER, site=['4fold']),
        expand(rules.angsd_estimate_sfs_byCity_byHabitat_bootstrap.output, habitat=HABITATS, city=CITIES, site=['4fold']),
        expand(rules.split_bootstrapped_sfs_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, number=NUMBER, site=['4fold']),
        expand(rules.angsd_bootstrapped_fst_readable.output, city=CITIES, number=NUMBER, site=['4fold']),
        expand(rules.angsd_bootstrapped_diversity_neutrality_stats_byCity_byHabitat.output, habitat=HABITATS, city=CITIES, number=NUMBER, site=['4fold'])  
    output:
        '{0}/angsd_byCity_byHabitat_bootstrap.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """

