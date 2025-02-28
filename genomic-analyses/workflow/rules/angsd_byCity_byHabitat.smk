# Rules for estimating SFS (1D and 2D), summary stats, and GLs for urban and rural habitats within cities

###############################
#### SFS AND SUMMARY STATS ####
###############################

rule angsd_saf_likelihood_byCity_byHabitat:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each habitat in each city using ANGSD. 
    Uses only 4fold sites.
    """
    input:
        unpack(get_files_for_saf_estimation_byCity_byHabitat)
    output:
        saf = '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.gz'.format(ANGSD_DIR),
        saf_idx = '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.idx'.format(ANGSD_DIR),
        saf_pos = '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.pos.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_saf_likelihood_byCity_byHabitat/{city}_{habitat}_{site}_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}'.format(ANGSD_DIR)
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '03:00:00'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*60/100 ));
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -baq 2 \
            -ref {input.ref} \
            -sites {input.sites} \
            -minInd $MIN_IND \
            -minMapQ 30 \
            -minQ 20 \
            -doSaf 1 \
            -anc {input.ref} \
            -bam {input.bams} 2> {log}
        """

rule angsd_estimate_joint_sfs_byCity:
    """
    Estimated folded, two-dimensional urban-rural SFS for each city using realSFS. Uses 4fold sites.
    """
    input:
        saf = get_habitat_saf_files_byCity,
        sites = rules.convert_sites_for_angsd.output,
        idx = rules.angsd_index_degenerate_sites.output,
    output:
        '{0}/sfs/by_city/{{city}}/{{city}}_{{site}}_r_u.2dsfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_2dsfs_byCity/{city}_{site}.2dsfs.log'
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
            -P {threads} > {output} 2> {log}
        """

rule angsd_estimate_sfs_byCity_byHabitat:
    """
    Estimate folded SFS separately for each habitat in each city (i.e., 1D SFS) using realSFS. 
    """
    input:
        saf = rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx,
        sites = rules.convert_sites_for_angsd.output,
        idx = rules.angsd_index_degenerate_sites.output,
    output:
        '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.sfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_sfs_byCity_byHabitat/{city}_{habitat}_{site}_sfs.log'
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
            -seed 42 > {output} 2> {log}
        """

#######################
#### FST AND THETA ####
#######################

rule angsd_fst_index:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Fst estimation.
    """
    input: 
        saf_idx = get_habitat_saf_files_byCity,
        joint_sfs = rules.angsd_estimate_joint_sfs_byCity.output,
        sites = rules.convert_sites_for_angsd.output
    output:
        fst = '{0}/summary_stats/hudson_fst/{{city}}/{{city}}_{{site}}_r_u.fst.gz'.format(ANGSD_DIR),
        idx = '{0}/summary_stats/hudson_fst/{{city}}/{{city}}_{{site}}_r_u.fst.idx'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_fst_index/{city}_{site}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    params:
        fstout = '{0}/summary_stats/hudson_fst/{{city}}/{{city}}_{{site}}_r_u'.format(ANGSD_DIR)
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

rule angsd_fst_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        rules.angsd_fst_index.output.idx
    output:
        '{0}/summary_stats/hudson_fst/{{city}}/{{city}}_{{site}}_r_u_readable.fst'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_fst_readable/{city}_{site}_readable_fst.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

rule angsd_estimate_thetas_byCity_byHabitat:
    """
    Generate per-site thetas in each habitat for each city from 1DSFS
    """
    input:
        saf_idx = rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_byCity_byHabitat.output,
        sites = rules.convert_sites_for_angsd.output
    output:
        idx = '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.thetas.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_thetas_byCity_byHabitat/{city}_{habitat}_{site}_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 4
    params:
        out = '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}'.format(ANGSD_DIR)
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

rule angsd_diversity_neutrality_stats_byCity_byHabitat:
    """
    Estimate pi, Waterson's theta, Tajima's D, etc. in each habitat in each city.
    """
    input:
        rules.angsd_estimate_thetas_byCity_byHabitat.output.idx
    output:
       '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_diversity_neutrality_stats_byCity_byHabitat/{city}_{habitat}_{site}_diversity_neutrality.log'
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

rule angsd_byCity_byHabitat_done:
    """
    Generate empty flag file signalling successful completion of SFS, summary stat and GL estimation 
    for habitats within cities
    """
    input:
        expand(rules.angsd_saf_likelihood_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, site=['4fold']),
        #expand(rules.angsd_estimate_joint_sfs_byCity.output, city=CITIES, site=['4fold']),
        #expand(rules.angsd_estimate_sfs_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, site=['4fold']),
        #expand(rules.angsd_fst_readable.output, city=CITIES, site=['4fold']),
        #expand(rules.angsd_diversity_neutrality_stats_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, site=['4fold'])
    output:
        '{0}/angsd_byCity_byHabitat.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """

